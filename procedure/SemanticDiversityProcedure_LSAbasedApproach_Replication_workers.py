# Import used modules
import os
import time
import random
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from sklearn.metrics.pairwise import cosine_similarity
from nltk.corpus.reader.bnc import BNCCorpusReader
from Context_LSA_SemanticDiversityProcedure_Replication_utils import TextPreProcessing, SplitDocuments, SemanticDiversity


def getContexts(lemmatization, path, wordWindow=1000):
    """ Returns a list of contexts (subdivided documents based on word window, default 1,000) of the British National Corpus.  

    Parameters
    ----------
    lemmatization : Boolean
        .
    path : str
        The directory of the folder where the BNC is saved.
    wordWindow : int
        The word length for which the BNC documents are subdivided.
        If not given, 1,000 is used as default. 

    Returns
    -------
    contextsList : list
        List of all the contexts in the corpus.
    contextInfo : str
        Id of the context.
    docInfo : list
        Id of the document.
    
    Notes
    -----

    References
    ----------

    """
    # Set time
    start_time = time.clock()
    # Set BNC reader & parameters
    bnc_reader = BNCCorpusReader(root="Resources/Corpora/BNC/Texts", fileids=r'[A-K]\/w*\/w*/.xml')
    # Check if text is from written source
    tags = [elem.tag for event, elem in ET.iterparse(path, events=("start", "end"))]
    if "wtext" in tags: 
        docID = path[33:-4]
        # Set time
        start_time = time.clock()
        # Read in a document as list of words
        docWordListRaw = bnc_reader.words(fileids=path[28:], strip_space=True, stem=lemmatization)
        # Preprocessing of raw text
        docWordList = TextPreProcessing(docWordListRaw)   
        # Split document into contexts 
        contextsList, contextInfo, docInfo = SplitDocuments(docWordList, docID, wordWindow) 
    elif "stext" in tags:
        contextsList = "SPOKEN"
        contextInfo = "SPOKEN"
        docInfo = "SPOKEN"
    else:
        contextsList = "NEITHER"
        contextInfo = "NEITHER"
        docInfo = "NEITHER"
    # Print out status
    t = time.clock()
    print('t: ', t/60, end='\t')
    print(t - start_time, "multiprocessor seconds")
    return (contextsList, contextInfo, docInfo)

def SaveDiversityInfo (contextVectors, wordInfo, maxDocsNumber=2000, svUnweighting=False):
    """ 
    Computes a word's semantic diversity. 

    Parameters
    ----------
    contextVectors : array-like
        The context vectors of all the 1,000-word contexts identified. 
        from which calculate a word's semantic 
        diversity.
    wordInfo: tuple
        wordInfo[0] : str
            The target word.
        wordInfo[1] : str
            The identification number of the target word.
        wordInfo[2] : list 
            A list of all the context labels in which the target word
            occurs in. 
    maxContextNumber : int, optional
        The maximum number of context vectors above which randomly 
        sample them in the calculation of the cosine similarity 
        matrix. If not give, 2,000 is used as maximum number of 
        vectors. 

    Returns
    -------

    Notes
    -----

    References
    ----------

    """
    
    # Set time
    start_time = time.clock()
    # Get word id and context list
    word, wordID, wordContextIndexes = wordInfo[0], wordInfo[1], wordInfo[2]
    wordContextIndexes = list(wordContextIndexes)
    # Randomly select maximum number of contexts if more than maximum given
    if len(wordContextIndexes)>maxDocsNumber:
        wordContextIndexes = random.sample(wordContextIndexes,maxDocsNumber)
    # Subset all of the context vectors in which the target word occurs only
    wordContextVectors = contextVectors[wordContextIndexes,:]
    # Compute semantic diversity of the target word
    if svUnweighting==True:
        method=="Hoffman"
    else:
        method=="hsiao"
    meanCosines, semanticDiversity = SemanticDiversity(word, wordContextVectors, method)
    # Print out status
    t = time.clock()
    print(word)
    print('t: ', t/60, end='\t')
    print(t - start_time, "multiprocessor seconds")
    return {'Word': word, 'SemanticDiversity': semanticDiversity ,'MeanCosines': meanCosines}
