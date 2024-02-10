from Bio import Entrez
from Bio import Medline
import functions_framework
import os
import openai
import arxiv
import time
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

@functions_framework.http

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

      # pprint.pprint(records)

      try:
          getPapers.append({
              "Title":records[0]["TI"],
              "Authors":records[0]["AU"],
              "Source":records[0]["SO"].split(". doi: ")[0],
              "Abstract":records[0]["AB"],
              "PMID":records[0]["PMID"],
              "Affiliation":records[0]["AD"],
              "DOI":records[0]["SO"].split(". doi: ")[1].split(".")[0],
              "paperLink":"https://pubmed.ncbi.nlm.nih.gov/" + records[0]["PMID"] + "/"
                            })
      except KeyError:
          getPapers.append({
              "Title":records[0]["TI"],
              "Authors":records[0]["AU"],
              "Source":records[0]["SO"].split(". doi: ")[0],
              "Abstract":"",
              "PMID":records[0]["PMID"],
              "Affiliation":records[0]["AD"],
              "DOI":records[0]["SO"].split(". doi: ")[1].split(".")[0],
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

    affiliation_list = [""] * len(authorsTmp)

    getPapers.append({
      "Title":result.title,
      "Authors":authorsTmp,
      "Source":"arXiv",
      "Abstract":result.summary,
      "PMID":"",
      "Affiliation":affiliation_list,
      "DOI":result.doi,
      "paperLink":result.links[0]
    })

  return(getPapers)

def bioRxiv_fetch_data(url):
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()['collection']
    else:
        return []

def contains_keywords(paper, keywords):
    title = paper['title'].lower()
    abstract = paper['abstract'].lower()
    return any(keyword in title or keyword in abstract for keyword in keywords)


def bioRxiv_multiPaperInfoGet(getPapers, termStr, bioArxivReldateStr):

  # Define the number of days back you want to fetch papers from
  # N_days_back = 1  # Example: Set to 1 for papers from 1 day ago
  N_days_back = int(bioArxivReldateStr)  # Example: Set to 1 for papers from 1 day ago

  # Calculate start and end dates
  end_date = datetime.datetime.now() - datetime.timedelta(days=N_days_back)
  start_date = end_date - datetime.timedelta(days=1)

  # Format dates in YYYY-MM-DD format for the URL
  start_date_str = start_date.strftime('%Y-%m-%d')
  end_date_str = end_date.strftime('%Y-%m-%d')

  # URLs for the API calls
  url1 = f"https://api.biorxiv.org/details/biorxiv/{start_date_str}/{end_date_str}/0/json"
  url2 = f"https://api.biorxiv.org/details/biorxiv/{start_date_str}/{start_date.strftime('%Y-%m-%d')}/0/json"

  # Define your keywords here
  keywords = termStr.lower().split(" ")
  print(keywords)

  # Fetch data from both URLs
  data_day1 = bioRxiv_fetch_data(url1)
  data_day2 = bioRxiv_fetch_data(url2)

  # Extract DOIs from the second day's data
  dois_day2 = {paper['doi'] for paper in data_day2}

  # Filter out papers from the first day that are not in the second day's data
  unique_papers = [paper for paper in data_day1 if paper['doi'] not in dois_day2]

  # Apply keyword search to the unique papers
  keyword_filtered_papers = [paper for paper in unique_papers if contains_keywords(paper, keywords)]

  for result in keyword_filtered_papers:
    authorsTmp=[]
    affiliation_list=[]
    authorsTmp = result["authors"].split("; ")
    affiliation_list = [""] * (len(authorsTmp) - 1) + [result["author_corresponding_institution"]]

    getPapers.append({
      "Title":result["title"],
      "Authors":authorsTmp,
      "Source":"bioRxiv",
      "Abstract":result["abstract"],
      "PMID":"",
      "Affiliation":affiliation_list,
      "DOI":result["doi"],
      "paperLink": "https://doi.org/" + result["doi"]
  })
  
  return(getPapers)

def chemRxiv_fetch_data(url):
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()['itemHits']
    else:
        return []

def chemRxiv_multiPaperInfoGet(getPapers, chemRxivtermStr, chemRxivReldateStr):
  # Define the number of days back you want to fetch papers from
  # N_days_back = 7  # Example: Set to 1 for papers from 1 day ago
  N_days_back = int(chemRxivReldateStr)

  # Calculate start and end dates
  end_date = datetime.datetime.now() - datetime.timedelta(days=N_days_back)
  start_date = end_date - datetime.timedelta(days=0)

  # Format dates in YYYY-MM-DD format for the URL
  start_date_str = start_date.strftime('%Y-%m-%d')
  end_date_str = end_date.strftime('%Y-%m-%d')

  # Set search strict
  search_strict=chemRxivtermStr

  # API URL
  url = f"https://chemrxiv.org/engage/chemrxiv/public-api/v1/items?term=%22{search_strict}%22&limit=50&searchDateFrom={start_date_str}&searchDateTo={end_date}"

  # Fetch data from the URL
  results = chemRxiv_fetch_data(url)

  # getPapers=[]
  for item in results:
    paper = item['item']

    authorNameList=[]
    institutionList=[]
    for author_names in paper['authors']:
      authorNameList.append(f"{author_names['firstName']} {author_names['lastName']}")
      if author_names['institutions']:
        institutionList.append(f"{author_names['institutions'][0]['name']}")
      else:
        institutionList.append(f"")

    getPapers.append({
      "Title":paper["title"],
      "Authors":authorNameList,
      "Source":"chemRxiv",
      "Abstract":paper["abstract"],
      "PMID":"",
      "Affiliation":institutionList,
      "DOI":paper["doi"],
      "paperLink": "https://doi.org/" + paper["doi"]
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
      time.sleep(0.5)
    else:
      getpaper["SummaryAbstract"] = ""

  return(getPapers)

def uploadDateTimeJST():
  return( datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))) )

def format_authors_institutions_with_missing(authors, institutions):
    formatted_authors = []  # To store the formatted author strings
    for author, inst in zip(authors, institutions):
        if inst:  # If institution information is present
            formatted_authors.append(f"{author} ({inst})")
        else:  # If institution information is missing
            formatted_authors.append(author)
    
    if len(authors) <= 3:
        return ', '.join(formatted_authors)
    else:
        # For more than 3 authors, include the first two and the last author in the output
        return f"{formatted_authors[0]}, {formatted_authors[1]} ... {formatted_authors[-1]}"

def uploadToTeams(getPapers, teamsWebHook):
  myTeamsMessage = pymsteams.connectorcard(teamsWebHook)
  for i, paperInfo in enumerate(getPapers):
    print("Upload_For_Teams:", i, paperInfo["Source"],paperInfo["Title"])
    AuthorsName = format_authors_institutions_with_missing(paperInfo["Authors"], paperInfo["Affiliation"])
    uploadDateTimeAndNum = uploadDateTimeJST().strftime('%Y/%m/%d %H:%M') + " | " + paperInfo["Source"] + " | " + "{0}/{1}".format(str(i+1), str(len(getPapers)))

    message = "**{0}** <br> **[{1}]({2})** <br> {3} <br><br> {4} <br><br> {5}".format(uploadDateTimeAndNum, paperInfo["Title"], paperInfo["paperLink"], AuthorsName, paperInfo["SummaryAbstract"], "")

    myTeamsMessage.text(message)
    myTeamsMessage.send()
  
def uploadToLINE(getPapers, LINE_token, LINE_channelID):
  line_bot_api = LineBotApi(LINE_token)
  for i, paperInfo in enumerate(getPapers):
    print("Upload_For_LINE:", i, paperInfo["Source"],paperInfo["Title"])
    AuthorsName = format_authors_institutions_with_missing(paperInfo["Authors"], paperInfo["Affiliation"])
    uploadDateTimeAndNum = uploadDateTimeJST().strftime('%Y/%m/%d %H:%M') + " | " + paperInfo["Source"] + " | " + "{0}/{1}".format(str(i+1), str(len(getPapers)))

    message = "{0}\n{1}\n{2}\n{3}\n\n{4}".format(uploadDateTimeAndNum, paperInfo["Title"], paperInfo["paperLink"], AuthorsName, paperInfo["SummaryAbstract"])

    try:
        line_bot_api.push_message(LINE_channelID, TextSendMessage(text=message, title=paperInfo["Title"]))
    except LineBotApiError as e:
        print("error")


def uploadToSlack(getPapers, SlackWebHookURL):
  for i, paperInfo in enumerate(getPapers):
    print("Upload_For_Slack:", i, paperInfo["Source"],paperInfo["Title"])
    AuthorsName = format_authors_institutions_with_missing(paperInfo["Authors"], paperInfo["Affiliation"])
    uploadDateTimeAndNum = uploadDateTimeJST().strftime('%Y/%m/%d %H:%M') + " | " + paperInfo["Source"] + " | " + "{0}/{1}".format(str(i+1), str(len(getPapers)))

    message = "*{0}*\n\n <{2}|*{1}*>\n {3} \n\n {4}".format(uploadDateTimeAndNum, paperInfo["Title"], paperInfo["paperLink"], AuthorsName, paperInfo["SummaryAbstract"])

    requests.post(SlackWebHookURL, data=json.dumps({
        "text" : message,
    }))


def uploadToDiscord(getPapers, DISCORD_WEBHOOK_URL):
  for i, paperInfo in enumerate(getPapers):
    print("Upload_For_Discord:", i, paperInfo["Source"],paperInfo["Title"])
    AuthorsName = format_authors_institutions_with_missing(paperInfo["Authors"], paperInfo["Affiliation"])
    uploadDateTimeAndNum = uploadDateTimeJST().strftime('%Y/%m/%d %H:%M') + " | " + paperInfo["Source"] + " | " + "{0}/{1}".format(str(i+1), str(len(getPapers)))

    message = "**{0}**\r**[{1}](<{2}>)**\r{3}\r\r{4}\r\r{5}".format(uploadDateTimeAndNum, paperInfo["Title"], paperInfo["paperLink"], AuthorsName, paperInfo["SummaryAbstract"], "​")

    discord = Discord(url=DISCORD_WEBHOOK_URL)
    discord.post(content = message,
                 username="LA_Article_Deliver")


def main(request):
  ENTREZ_EMAIL = keys_tokens.ENTREZ_EMAIL
  teamsWebHook = keys_tokens.TEAMS_WEBHOOK
  SlackWebHookURL = keys_tokens.SLACK_WEBHOOK_URL
  OPENAI_API_KEY = keys_tokens.OPENAI_API_KEY
  LINE_token = keys_tokens.LINE_TOKEN
  LINE_channelID = keys_tokens.LINE_CHANNEL_ID
  DISCORD_WEBHOOK_URL = keys_tokens.DISCORD_WEBHOOK_URL

  termStr = "<YOUR SEARCH WORD HERE>" # ex.) lab automation
  termStrForBioRxiv="<YOUR SEARCH WORD HERE>" # ex.)automation
  termStrForChemRxiv="<YOUR SEARCH WORD HERE>" # ex.)automation
  datetypeStr = "edat" #'mdat' (modification date), 'pdat' (publication date) and 'edat' (Entrez date)
  reldateStr = "1" # The search range will be set from today back to N days. The minimum value is 1.
  arxivReldateStr = "3" # The search range will be set from today back to N days. The database update is slow, so the minimum value is 3. 
  bioRxivReldateStr = "1" # The search range will be set from today back to N days. The minimum value is 1. 
  chemRxivReldateStr = "1" # The search range will be set from today back to N days. The minimum value is 1. 

  getPapers=[]
  # to get papers from pubmed
  pubmed_record = pubmed_searchToGetIDs(termStr, datetypeStr, reldateStr, ENTREZ_EMAIL)
  papers = pubmed_medline_multiPaperInfoGet(getPapers, pubmed_record, ENTREZ_EMAIL)
  print("Pubmed papers: ",len(papers))

  # to get papers from bioRxiv
  papers = bioRxiv_multiPaperInfoGet(papers, termStrForBioRxiv, bioRxivReldateStr)
  print("Pubmed + bioRxiv papers: ",len(papers))

  # to get papers from chemRxiv
  papers = chemRxiv_multiPaperInfoGet(papers, termStrForChemRxiv, chemRxivReldateStr)
  print("Pubmed + bioRxiv + chemRxiv papers: ",len(papers))

  # to get papers from arxiv
  arxiv_termStr=termStr
  papers = arxiv_multiPaperInfoGet(papers, arxiv_termStr, arxivReldateStr)
  print("Pubmed + bioRxiv + chemRxiv + arXiv papers (Total) : ",len(papers))

  # summarize the abstract using GPT model
  #basePrompt = "以下は英語学術論文のアブストラクトです。本論文の主張を、端的に日本語で、100-200文字で示していただけますか? なお、日本語的に不可解な単語は無理に訳さずに英語のまま表現してください。"
  basePrompt = "以下は英語学術論文のアブストラクトです。本論文の主張を、端的に示す10個程度の日本語の単語を\", \"で区切って示していただけますか？なお、日本語的に不可解な単語は無理に訳さずに英語の単語のまま表現してください。"
  # model = "gpt-3.5-turbo"
  # model = "gpt-4-1106-preview"	
  model = "gpt-4-0125-preview"	# Select favorid OpenAI model.
  papers = abstractGPTsummarize(papers, OPENAI_API_KEY, basePrompt, model)

  # upload the paper info to Teams
  uploadToTeams(papers, teamsWebHook)

  # send notification to LINE
  uploadToLINE(papers, LINE_token, LINE_channelID)

  # upload the paper info to Slack
  uploadToSlack(papers, SlackWebHookURL)

  # upload the paper info to Discord
  uploadToDiscord(papers, DISCORD_WEBHOOK_URL)

  return "All Done."


if __name__ == "__main__":
  main()
