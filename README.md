
![Overall](????????)

## English
# Automated-Paper-Summary-Bot

This Python script uses the Bio.Entrez library, OpenAI's GPT-3.5-turbo model, and the ArXiv API to automatically search and retrieve academic papers from PubMed and ArXiv, generate summaries of their abstracts, and share these summaries along with paper details to Microsoft Teams, Slack, and LINE.


## Features
- Automatically searches for academic papers from PubMed and ArXiv based on specific search terms.
- Summarizes these papers' abstracts using the OpenAI GPT-3.5-turbo model (default).
- Shares the summarized abstracts along with paper details (title, authors, source, abstract, link, etc.) on Microsoft Teams, Slack, and LINE.
- If any paper source (Pubmed, ArXiv) or sharing destination (Teams, Slack, LINE) is not necessary, please comment them out to deal with them.

## How to Use
1. Download or clone this repository.
2. Update the main() function with the desired search term, date type, and relative date.
3. Make sure to import your sensitive data (API keys, tokens, webhooks) from a separate file to avoid exposure.
4. Run the script by typing `python main.py` in your terminal/command prompt.

## Dependencies
- BioPython
- openai
- arxiv
- pymsteams
- line-bot-sdk
- requests

## Japanese
# 自動論文要約ボット
このPythonスクリプトは、Bio.Entrezライブラリ、OpenAIのGPT-3.5-turboモデル(など)、そしてArXiv APIを利用して、PubMedとArXivから自動的に学術論文を検索・取得し、それらのアブストラクトの要約を生成し、その要約と論文の詳細をMicrosoft Teams、Slack、そしてLINEに共有します。

## 機能
- 特定の検索語に基づいてPubMedとArXivから学術論文を自動検索します。
- これらの論文のアブストラクトをOpenAIのGPT-3.5-turboモデル(既定)を使って要約します。
- 要約したアブストラクトと共に論文の詳細（タイトル、著者、ソース、アブストラクト、リンクなど）をMicrosoft Teams、Slack、そしてLINEに共有します。
- 論文ソース元(Pubmed, ArXiv)および共有先(Teams, Slack, LINE)で不要なものはコメントアウトして対応してください。


## 使い方
- このリポジトリをダウンロードまたはクローンします。
- `main()`関数を、希望する検索語、日付タイプ、相対日付に更新します。
- 漏洩を避けるため、センシティブなデータ（APIキー、トークン、ウェブフック）を別のファイル(`keys_tokens.py`)からインポートするようにします。Token情報などをそちらに入力してください。
- ターミナル/コマンドプロンプトで`python main.py`と入力してスクリプトを実行します。


## 依存関係
- BioPython
- openai
- arxiv
- pymsteams
- line-bot-sdk
- requests