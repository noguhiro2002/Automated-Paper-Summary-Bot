from Bio import Entrez
from Bio import Medline
import os
import openai
import arxiv
import datetime
import pymsteams
import requests
import json
from linebot import LineBotApi
from linebot.models import TextSendMessage
from linebot.exceptions import LineBotApiError
from discordwebhook import Discord

# Enter your keys and tokens in a separate file and import it
# This can prevent accidental exposure of sensitive information
import keys_tokens

def pubmed_searchToGetIDs(termStr, datetypeStr, reldateStr, ENTREZ_EMAIL):
  Entrez.email = ENTREZ_EMAIL
  handle = Entrez.esearch(db="pubmed", term=termStr,  datetype=datetypeStr, reldate=reldateStr)
  record = Entrez.read(handle)
  return(record)

def pubmed_medline_multiPaperInfoGet(getPapers,pubmedRecord, ENTREZ_EMAIL):
  Entrez.email = ENTREZ_EMAIL
  record=pubmedRecord
  # getPapers=[]
  records=[]
  for ID in record["IdList"]:
      handle = Entrez.efetch(db="pubmed", id=ID, rettype="medline", retmode="text")
      records = Medline.parse(handle)
      records = list(records)

      try:
          getPapers.append({
              "Title":records[0]["TI"],
              "Authors":records[0]["AU"],
              "Source":records[0]["SO"],
              "Abstract":records[0]["AB"],
              "PMID":records[0]["PMID"],
              "DOI":"",
              "paperLink":"https://pubmed.ncbi.nlm.nih.gov/" + records[0]["PMID"] + "/"
                            })
      except KeyError:
          getPapers.append({
              "Title":records[0]["TI"],
              "Authors":records[0]["AU"],
              "Source":records[0]["SO"],
              "Abstract":"",
              "PMID":records[0]["PMID"],
              "DOI":"",
              "paperLink":"https://pubmed.ncbi.nlm.nih.gov/" + records[0]["PMID"] + "/"
                            })
  return(getPapers)

def arxiv_multiPaperInfoGet(getPapers, termStr, arxivReldateStr):
  search_minus_Days=int(arxivReldateStr)

  current_date = datetime.datetime.now()

  minusDaysForSearch=[]
  minusDaysForSearch.append((current_date - datetime.timedelta(days=search_minus_Days)).strftime("%Y%m%d"))
  minusDaysForSearch.append((current_date - datetime.timedelta(days=search_minus_Days-1)).strftime("%Y%m%d"))

  search = arxiv.Search(
    query = termStr + " AND submittedDate:[{0} TO {1}]".format(minusDaysForSearch[0], minusDaysForSearch[1]),
    id_list = [],
    max_results = 10,
    sort_by = arxiv.SortCriterion.Relevance,
    sort_order = arxiv.SortOrder.Descending
  )

  for result in search.results():
    authorsTmp=[]
    for author in result.authors:
      authorsTmp.append(str(author))

    getPapers.append({
      "Title":result.title,
      "Authors":authorsTmp,
      "Source":"arXiv",
      "Abstract":result.summary,
      "PMID":"",
      "DOI":result.doi,
      "paperLink":result.links[0]
    })

  return(getPapers)


def abstractGPTsummarize(getPapers, OPENAI_API_KEY, basePrompt, model):
  openai.api_key = OPENAI_API_KEY

  for i, getpaper in enumerate(getPapers):
    if getpaper["Abstract"] != "":
      Abstract = getpaper["Abstract"]
      messageForGPT = basePrompt + "　" + Abstract
      renponse=""
      renponse = openai.ChatCompletion.create(
        model=model,
        messages=[
            {"role":"system", "content":"You are a helpful assistant."},
            {"role": "user", "content": messageForGPT}
        ],
        temperature=1,
        max_tokens=1001,
        top_p=1,
        frequency_penalty=0,
        presence_penalty=0
      )
      paperGPTSumamry = renponse.choices[0].message.content.strip()
      getpaper["SummaryAbstract"] = paperGPTSumamry
    else:
      getpaper["SummaryAbstract"] = ""

  return(getPapers)

def uploadDateTimeJST():
  return( datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))) )

def uploadToTeams(getPapers, teamsWebHook):
  myTeamsMessage = pymsteams.connectorcard(teamsWebHook)

  for i, paperInfo in enumerate(getPapers):
    AuthorsName=""
    mainBody=""

    # Create upload AuthorsName
    for n, authorName in enumerate(paperInfo["Authors"]):
      if (n+1) != len(paperInfo["Authors"]):
        AuthorsName = AuthorsName + authorName + ", "
      else:
        AuthorsName = AuthorsName + authorName

    uploadDateTimeAndNum = uploadDateTimeJST().strftime('%Y/%m/%d %H:%M') + " | " + "{0}/{1}".format(str(i+1), str(len(getPapers)))

    mainBody = "**{0}** <br><br> **[{1}]({2})** <br> {3} <br> *{4}* <br> <br> **要約:** <br> {5} <br> <br> **Abstract:** <br> {6} <br> <br>".format(uploadDateTimeAndNum, paperInfo["Title"], paperInfo["paperLink"], AuthorsName, paperInfo["Source"], paperInfo["SummaryAbstract"], paperInfo["Abstract"])
    myTeamsMessage.text(mainBody)
    myTeamsMessage.send()
  
def uploadToLINE(getPapers, LINE_token, LINE_channelID):
  line_bot_api = LineBotApi(LINE_token)

  for i, paperInfo in enumerate(getPapers):
    AuthorsName=""
    mainBody=""
    for n, authorName in enumerate(paperInfo["Authors"]):
      if (n+1) != len(paperInfo["Authors"]):
        AuthorsName = AuthorsName + authorName + ", "
      else:
        AuthorsName = AuthorsName + authorName

    uploadDateTimeAndNum = uploadDateTimeJST().strftime('%Y/%m/%d %H:%M') + " | " + "{0}/{1}".format(str(i+1), str(len(getPapers)))

    mainBoey = "{0}\n\n{1}\n{2}\n{3}\n{4}\n{5}".format(uploadDateTimeAndNum, paperInfo["Title"], paperInfo["paperLink"], AuthorsName, paperInfo["Source"], paperInfo["SummaryAbstract"])
    try:
        line_bot_api.push_message(LINE_channelID, TextSendMessage(text=mainBoey, title=paperInfo["Title"]))
    except LineBotApiError as e:
        print("error")


def uploadToSlack(getPapers, SlackWebHookURL):
  for i, paperInfo in enumerate(getPapers):
    AuthorsName=""
    message=""
    for n, authorName in enumerate(paperInfo["Authors"]):
      if (n+1) != len(paperInfo["Authors"]):
        AuthorsName = AuthorsName + authorName + ", "
      else:
        AuthorsName = AuthorsName + authorName

    uploadDateTimeAndNum = uploadDateTimeJST().strftime('%Y/%m/%d %H:%M') + " | " + "{0}/{1}".format(str(i+1), str(len(getPapers)))
    message = "*{0}*\n\n <{2}|*{1}*>\n {3}\n _{4}_\n\n *要約:*\n {5}\n\n *Abstract:*\n{6}".format(uploadDateTimeAndNum, paperInfo["Title"], paperInfo["paperLink"], AuthorsName, paperInfo["Source"], paperInfo["SummaryAbstract"], paperInfo["Abstract"])

    requests.post(SlackWebHookURL, data=json.dumps({
        "text" : message,
    }))

def uploadToDiscord(getPapers, DISCORD_WEBHOOK_URL):
  for i, paperInfo in enumerate(getPapers):
    AuthorsName=""
    message=""
    for n, authorName in enumerate(paperInfo["Authors"]):
      if (n+1) != len(paperInfo["Authors"]):
        AuthorsName = AuthorsName + authorName + ", "
      else:
        AuthorsName = AuthorsName + authorName

    uploadDateTimeAndNum = uploadDateTimeJST().strftime('%Y/%m/%d %H:%M') + " | " + "{0}/{1}".format(str(i+1), str(len(getPapers)))
    # message = "**{0}**\r\r**{1}**\r{2}\r{3}\r_{4}_\r\r**要約:**\r{5}\r\r**Abstract:**\r{6}".format(uploadDateTimeAndNum, paperInfo["Title"], paperInfo["paperLink"], AuthorsName, paperInfo["Source"], paperInfo["SummaryAbstract"], paperInfo["Abstract"])
    message = "**{0}**\r\r**{1}**\r{2}\r{3}\r_{4}_\r\r要約:\r{5}".format(uploadDateTimeAndNum, paperInfo["Title"], paperInfo["paperLink"], AuthorsName, paperInfo["Source"], paperInfo["SummaryAbstract"])

    discord = Discord(url=DISCORD_WEBHOOK_URL)
    discord.post(content = message,
                 username="Postman") #you can change your username to post.

def main(data, context):
  ENTREZ_EMAIL = keys_tokens.ENTREZ_EMAIL
  teamsWebHook = keys_tokens.TEAMS_WEBHOOK
  SlackWebHookURL = keys_tokens.SLACK_WEBHOOK_URL
  OPENAI_API_KEY = keys_tokens.OPENAI_API_KEY
  LINE_token = keys_tokens.LINE_TOKEN
  LINE_channelID = keys_tokens.LINE_CHANNEL_ID
  DISCORD_WEBHOOK_URL = keys_tokens.DISCORD_WEBHOOK_URL

  termStr = "lab automation" # type your search term. ex) lab automation
  datetypeStr = "edat" #'mdat' (modification date), 'pdat' (publication date) and 'edat' (Entrez date)
  reldateStr = "1" # The search range will be set from today back to N days. The minimum value is 1.
  arxivReldateStr = "3" # The search range will be set from today back to N days. The database update is slow, so the minimum value is 3. 

  getPapers=[]
  # to get papers from pubmed
  pubmed_record = pubmed_searchToGetIDs(termStr, datetypeStr, reldateStr, ENTREZ_EMAIL)
  papers = pubmed_medline_multiPaperInfoGet(getPapers, pubmed_record, ENTREZ_EMAIL)

  # to get papers from arxiv
  arxiv_termStr=termStr
  papers = arxiv_multiPaperInfoGet(papers, arxiv_termStr, arxivReldateStr)

  # summarize the abstract using GPT model
  basePrompt = "以下は英語学術論文のアブストラクトです。本論文の主張を、端的に日本語で、100-200文字で示していただけますか? なお、日本語的に不可解な単語は無理に訳さずに英語のまま表現してください。" #type your prompt. Of course, you can change your language.
  model = "gpt-3.5-turbo" # select your OpenAI model
  papers = abstractGPTsummarize(papers, OPENAI_API_KEY, basePrompt, model)

  # # upload the paper info to Teams
  uploadToTeams(papers, teamsWebHook)

  # # send notification to LINE
  uploadToLINE(papers, LINE_token, LINE_channelID)

  # # upload the paper info to Slack
  uploadToSlack(papers, SlackWebHookURL)

  # upload the paper info to Discord
  uploadToDiscord(papers, DISCORD_WEBHOOK_URL)


if __name__ == "__main__":
  main()