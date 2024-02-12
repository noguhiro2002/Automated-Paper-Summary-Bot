
![Overall](https://github.com/noguhiro2002/Automated-Paper-Summary-Bot/blob/main/readme/Picture_1.png?raw=true)

## English
# Automated Paper Search/Summarize Bot
This Python script utilizes the Bio.Entrez library, OpenAI's GPT-4-turbo model, and the ArXiv API to automatically search and retrieve academic papers from PubMed, BioRxiv, ChemRxiv, and ArXiv. It generates about ten representative keywords from their abstracts and shares these keywords along with the paper details (title with link, authors, source, link, etc.) to Microsoft Teams, Slack, LINE, and Discord.

## Features
- Automatically searches for academic papers from PubMed, BioRxiv, ChemRxiv, and ArXiv based on specific search terms.
- Extracts about ten keywords that symbolize the paper from its abstract using the OpenAI GPT-4-turbo model (default).
- Shares the extracted keywords along with paper details (title (link), authors, source, link, etc.) on Microsoft Teams, Slack, LINE, and Discord.
- Feel free to use any search keywords you like.
- If any paper source (PubMed, BioRxiv, ChemRxiv, ArXiv) or sharing destination (Teams, Slack, LINE, Discord) is not necessary, please comment them out to deal with them.

## How to Use
1. Download or clone this repository.
2. Update the main() function with the desired search term, date type, and relative date.
3. Make sure to import your sensitive data (API keys, tokens, webhooks) from a separate file (keys_tokens.py) to avoid exposure. Enter your token information there.
4. Run the script by typing python main.py in your terminal/command prompt.
5. main_gcp.py is a script that operates with Google Cloud Platform (GCP)'s Cloud Function. Please feel free to use it.

## Dependencies
- BioPython
- openai
- arxiv
- pymsteams
- line-bot-sdk
- requests
- discordwebhook



## Japanese
# 自動論文検索/要約ボット
このPythonスクリプトは、Bio.Entrezライブラリ、OpenAIのGPT-4-turboモデル(など)、そしてPubMed/BioRxiv/ChemRxiv/ArXivのAPIを用いて、自動的に学術論文を検索・取得し、それらのアブストラクトから論文を象徴する単語を10個ほど生成し、それと論文の詳細をMicrosoft Teams、Slack、LINE、そしてDiscordに共有します。

## 機能
- 特定の検索語に基づいてPubMed/BioRxiv/ChemRxiv/ArXivから、特定のキーワードに関連する学術論文を自動検索します。
- これらの論文のアブストラクトから、その論文を象徴するキーワード１０個ほどをOpenAIのGPT-4-turboモデル(既定)を使って抽出します。
- 抽出したキーワードと共に論文の詳細（タイトル(リンク)、著者、ソース、リンクなど）をMicrosoft Teams、Slack、LINE、そしてDiscordに共有します。
- 検索キーワードはお好きなものをお使いください。
- 論文ソース元(Pubmed, BioRxiv, ChemRxiv, ArXiv)および共有先(Teams, Slack, LINE, Discord)で不要なものはコメントアウトして対応してください。

## 使い方
- このリポジトリをダウンロードまたはクローンします。
- `main()`関数を、希望する検索語、日付タイプ、相対日付に更新します。
- 漏洩を避けるため、センシティブなデータ（APIキー、トークン、ウェブフック）を別のファイル(`keys_tokens.py`)からインポートするようにします。Token情報などをそちらに入力してください。
- ターミナル/コマンドプロンプトで`python main.py`と入力してスクリプトを実行します。
- main_gcp.py は、Google Cloud Platform (GCP)のCloud Functionで動作するスクリプトです。是非お使いください。

## 依存関係
- BioPython
- openai
- arxiv
- pymsteams
- line-bot-sdk
- requests
- discordwebhook
