import requests
import json
import os

from ccm_demo.structure import structure
from ccm_demo.sequence import sequence
from ccm_demo.simulation import simulation
from ccm_demo.literature import literature

class Protein:
    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id
        self.json=self.gather()
        self.process_uniprot()
        self.description=self.get_description()
        self.references=self.extract_references()
        self.features=self.extract_features()
        self.comments=self.extract_comments()
        self.cross_references=self.get_xrefs()
        

    def gather(self):
        url="https://rest.uniprot.org/uniprotkb/{}?format=json".format(self.uniprot_id)
        response=requests.get(url)
        response.raise_for_status()
        
        content=response.content.decode().strip()
        content=json.loads(content)
        return content

    def process_uniprot(self):
        self.sequence=self.json["sequence"]["value"]
        self.organism={"name":self.json["organism"]["scientificName"],
                       "taxid":self.json["organism"]["taxonId"]}
        self.name=self.json["proteinDescription"]["recommendedName"]["fullName"]["value"]
        self.genes=[item["geneName"]["value"] for item in self.json["genes"]]
        self.comments=self.json["extraAttributes"]["countByCommentType"]
        self.features=self.json["extraAttributes"]["countByFeatureType"]
        
                       
    def extract_features(self, feature_types=None):
        if feature_types is not None:
            features=[feat for feat in self.json["features"] if feat["type"] in feature_types]
        else:
            features=[feat for feat in self.json["features"]]

        self.features=features
                      
        
    def extract_comments(self, comment_types):
        if comment_types is not None:
            comments=[feat for feat in self.json["comments"] if feat["type"] in comment_types]
        else:
            comments=[feat for feat in self.json["comments"]]

        self.comments=comments
        
    def extract_references(self):
        refs=[]
        for reference in self.json["references"]:
            ref={"title":reference["citation"]["title"],
                 "dbs":[{"database":item["database"], 
                         "id":item["id"]} for item in reference["citationCrossReferences"]]}
            refs.append(ref)
        self.references=refs

    def get_xrefs(self):
        xrefs=self.json["uniProtKBCrossReferences"]
        self.xref_types=list(set([item["database"] for item in xrefs]))
        self.xrefs=xrefs

    def get_cross_references(self, type=None):
        if type is None:
            return self.xrefs
        else:
            return [item for item in self.xrefs if item["database"]==type]


    def get_description(self):
        #TODO add summarize
        desc=[]
        for comment in self.json["comments"]:
            if "texts" in comment.keys():
                desc.append("\n".join([item["value"] for item in comment["texts"]]))

        self.description="\n".join(desc)
    
    def __str__(self):
        return self.description

    def __repr__(self):
         return "Protein instance with name {} and {} aa long".format(self.name, len(self.sequence))





