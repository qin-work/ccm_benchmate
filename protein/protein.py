import requests
import json
import pandas as pd
import os


class Protein:
    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id
        self.json=self.gather(self.uniprot_id)
        self.process_uniprot()
        

    def gather(self):
        url="https://rest.uniprot.org/uniprotkb/{}?format=json".format(self.uniprot_id)
        response=requests.get(url)
        response.raise_for_status()
        
        content=response.content.decode().strip()
        content=json.loads(content)
        return content

    def process_uniprot(self):
        self.sequence=self.json["sequence"]["value"]
        self.organism={"name":self.json["organism"]["scientififName"], 
                       "taxid":self.json["organism"]["taxonId"]}
        self.name=self.json["proteinDescription"]["recommendedName"]["fullName"]["value"]
        self.genes=[item["geneName"]["value"] for item in self.json["genes"]]
        self.comments=self.json["extraAttributes"]["countcountByCommentType"]
        self.features=self.json["extraAttributes"]["countcountByFeatureType"]
        
                       
    def extract_features(self, feature_types=None):
        if feature_types is not None:
            features=[feat for feat in self.json["features"] if feat["type"] in feature_types]
        else:
            features=[feat for feat in self.json["features"]]

        self.features=features
        return self
                      
        
    def extract_comments(self, comment_types):
        if comment_types is not None:
            comments=[feat for feat in self.json["comments"] if feat["type"] in comment_types]
        else:
            comments=[feat for feat in self.json["comments"]]

        self.features=comments
        return self
        
    def extract_references(self):
        refs=[]
        for reference in self.json["references"]:
            ref={"title":reference["title"], 
                 "dbs":[{"database":item["database"], 
                         "id":item["id"]} for item in reference["citationCrossReferences"]]}
            refs.append(ref)
        self.references=refs
        return self

    def get_description(self):
        #TODO add summarize
        desc=[]
        for comment in self.json["comments"]:
            if["texts"] in comment.keys():
                desc.append("\n".join([item["value"] for item in comment["texts"]]))

        self.description="\n".join(desc)
        return self
    
    def __str__(self):
        if self.description is not None:
            print(self.description)
        else: 
            print(self.get_description())

    def __repr__(self):
        print("Protein instance with name {} and {} aa long".format(self.name, len(self.sequence)))

        



