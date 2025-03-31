from bs4 import BeautifulSoup
import requests


class LiteratureSearch:
    def __init__(self, pubmed_url, arxiv_url, pubmed_api_key):
        self.pubmed_url = pubmed_url
        self.arxiv_url = arxiv_url
        self.pubmed_key = pubmed_api_key

    def search(self, query, database="pubmed"):
        pass



class Paper:
    def __init__(self, id):
        self.id = id
        self.title = None
        self.abstract = None
        self.text = None
        self.figures = None
        self.tables = None
        self.references = None

    def get_meta(self, source, get_references=False, get_keywords=False):
        if source =="pubmed":
            response=requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={}".format(self.id))
            response.raise_for_status()
            soup=BeautifulSoup(response.text, "xml")
            self.title = soup.find("ArticleTitle").text
            self.abstract = soup.find("AbstractText").text
        elif source == "arxiv":
            response = requests.get("http://export.arxiv.org/api/query?search_query=id:{}".format(self.id))
            response.raise_for_status()
            soup=BeautifulSoup(response, "xml")
            self.title = soup.find("title").text
            self.abstract = soup.find("summary").text
        else:
            raise NotImplementedError("source must be pubmed or arxiv other sources are not implemented")

    def download(self):
        #TODO get html/pdf of the paper ideally former
        pass

    def process(self, pdf):
        pass


class KnowledgeBase:
    def __init__(self, name, papers):
        self.name = name
        self.papers = papers
        self.text_embedding_model=None
        self.image_embedding_model=None

    def build(self, format="pkl"):
        pass

    def RAG(self):
        pass



