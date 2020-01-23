# Import used modules
import os
import re
import glob
import random
import numpy as np
import pandas as pd
from nltk.corpus import stopwords
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import cosine_similarity
import scipy.sparse as sp
from sklearn.preprocessing import normalize
from sklearn.utils.validation import check_is_fitted


def TextPreProcessing(rawWordList):
    """  """
    # Convert to lower case
    wordListLower = [w.lower() for w in rawWordList]
    # Non-alphabetic characters splitting
    wordListClean = [subword for word in wordListLower for subword in re.split('[^a-zA-Z]',word) if subword!=""]
    # Remove stop words
    stop_words = stopwords.words('english')
    wordListClean = [word for word in wordListClean if word not in stop_words]
    # Remove words with length <= 1
    wordListClean = [word for word in wordListClean if len(word)>1]
    return (wordListClean)

def SplitDocuments(docWordList, docID, contextWindow):
    """  """
    contextList = []
    contextInfo = []
    docInfo = []
    # Number of contexts in the document 
    nWords = len(docWordList)
    nContexts = nWords//contextWindow
    # Loop through each context and save dictionary of counts
    for nContext in range(0,nContexts):
        # Select context 
        initialIndex = contextWindow*nContext
        finalIndex = contextWindow*(1+nContext)
        contextWordList = docWordList[initialIndex:finalIndex]
        contextID = docID+"_"+str(nContext+1)
        contextInfo.append(contextID)
        docInfo.append(docID)
        # Append to list of all the contexts 
        contextList.append(" ".join(contextWordList))  
        
    return(contextList, contextInfo, docInfo)

def SemanticDiversity (word, wordContextVectors):
    """ 
    Returns a word's semantic diversity as the mean semantic distance
    between the context vectors in which the word occurs. 
    
    Semantic diversity can be calculated using either of the following
    methods, specified in the `method` parameter. The 'hoffman' method 
    uses the negative mean log to the base 10 cosine distances between 
    the context vectors described in [1], whereas the 'hsiao' method 
    uses the negative mean log to the base 2 cosine distances between 
    the context vectors as described in [2].

    Parameters
    ----------
    word : str
        The word for which to calculate semantic diversity.
    contextVectors : array-like
        The context vectors from which calculate a word's semantic 
        diversity.     
    method : {'hoffman', 'hsiao'}, 
        The method to use for the calculation of semantic diversity. 
    

    Returns
    -------
    meanCosines: float
        The mean of the pairwise cosine distances of all of the context 
        vectors in which the word occurs. 
    semanticDiversity: float
        The negative log of the mean cosine distances, aka a word's 
        semantic diversity.

    Notes
    -----


    References
    ----------
    .. [1] Hoffman, P., Lambon Ralph, M. A., & Rogers, T. T. (2013).
    Semantic diversity: A measure of semantic ambiguity based on
    variability in the contextual usage of words. *Behavior Research
    Methods, 45(3)*, 718730. doi:10.3758/s13428-012-0278-x
    .. [2] Hsiao, Y. & Nation, K. (2018). Semantic Diversity, 
    Frequency and the Development of Lexical Quality in Children’s 
    Word. *Journal of Memory and Language, 103*, 114-126.

    """
       
    # Compute pairwise cosine similarity between all of the context vectors
    cosineSimilarityMatrix = cosine_similarity(wordContextVectors, dense_output=False)
    # Save pairwise similarity matrix
    #np.savetxt(dataFolder+"PairwiseSimilarities\\"+word+".csv", cosSimilarities, delimiter=",")
    #Select the upper triangle of the similarity matrix only
    upperTriangle = cosineSimilarityMatrix[np.triu_indices(cosineSimilarityMatrix.shape[0],1)]
    # Compute mean of the upper triangle of the cosine similarity matrix
    meanCosines = upperTriangle.mean()
    # Compute a word's semantic diversity    
    if method=="hoffman":
        semanticDiversity = -np.log10(meanCosines)
    elif method=="hsiao":
        semanticDiversity = -np.log2(meanCosines)
    
    return (meanCosines, SemanticDiversity)


   
class LogEntropyVectorizer(CountVectorizer):
    # adatped from https://github.com/titipata/science_concierge/blob/master/science_concierge/vectorizer.py
    """Log-entropy vectorizer
    Convert collection of raw documents to matrix of log-entropy features
    Adds on functionality for scikit-learn CountVectorizer to
    calculate log-entropy term matrix
    Log-entropy
    -----------
    Assume we have term i in document j can be calculated as follows
    Global entropy
        p_ij = f_ij / sum_j(f_ij)
        g_i = 1 + sum_j (p_ij * log p_ij / log n)
    log-entropy of term i in document j is
        l_ij = log(1 + f_ij) * g_i
    where
        f_ij is number of term i that appears in document j
        sum_j(f_ij) is total number of times term i occurs in
            the whole documents
        n is total number of documents
        g_i is sum of entropy across all documents j
    Parameters
    ----------
    encoding : string, 'utf-8' by default.
        If bytes or files are given to analyze, this encoding is used to
        decode.
    decode_error : {'strict', 'ignore', 'replace'}
        Instruction on what to do if a byte sequence is given to analyze that
        contains characters not of the given `encoding`. By default, it is
        'strict', meaning that a UnicodeDecodeError will be raised. Other
        values are 'ignore' and 'replace'.
    ngram_range : tuple (min_n, max_n)
        The lower and upper boundary of the range of n-values for different
        n-grams to be extracted. All values of n such that min_n <= n <= max_n
        will be used.
    stop_words : string {'english'}, list, or None (default)
    lowercase : boolean, default True
        Convert all characters to lowercase before tokenizing.
    token_pattern : string
        Regular expression denoting what constitutes a "token", only used
        if ``analyzer == 'word'``. The default regexp selects tokens of 2
        or more alphanumeric characters (punctuation is completely ignored
        and always treated as a token separator).
    max_df : float in range [0, 1] or int, default=1.0
    min_df : float in range [0, 1] or int, default=1
    norm : 'l1', 'l2' or None, optional
        Norm used to normalize term vectors. None for no normalization.
    smooth_idf: boolean, default=False
    See also
    --------
    CountVectorizer
        Tokenize the documents and count the occurrences of token and return
        them as a sparse matrix
    TfidfTransformer
        Apply Term Frequency Inverse Document Frequency normalization to a
        sparse matrix of occurrence counts.
    Example
    -------
    >> model = LogEntropyVectorizer(norm=None, ngram_range=(1,1))
    >> docs = ['this this this book',
               'this cat good',
               'cat good shit']
    >> X = model.fit_transform(docs)
    References
    ----------
        - https://en.wikipedia.org/wiki/Latent_semantic_indexing
        - http://webpages.ursinus.edu/akontostathis/KontostathisHICSSFinal.pdf
    """
    def __init__(self, encoding='utf-8', decode_error='strict',
                 lowercase=True, preprocessor=None, tokenizer=None,
                 analyzer='word', stop_words=None, token_pattern=r"(?u)\b\w\w+\b",
                 vocabulary=None, binary=False,
                 ngram_range=(1, 1), max_df=1.0, min_df=1, min_tf=1,
                 max_features=None, norm='l2', smooth_idf=False):


        super(LogEntropyVectorizer, self).__init__(
            encoding=encoding,
            decode_error=decode_error,
            lowercase=lowercase,
            preprocessor=preprocessor,
            tokenizer=tokenizer,
            analyzer=analyzer,
            stop_words=stop_words,
            token_pattern=token_pattern,
            ngram_range=ngram_range,
            max_df=max_df,
            min_df=min_df,
            min_tf=min_tf,
            max_features=max_features,
            vocabulary=vocabulary,
            binary=binary,
        )

        self.norm = norm
        self.smooth_idf = smooth_idf


    def fit(self, raw_documents, y=None):
        """Learn vocabulary and log-entropy from training set.
        Parameters
        ----------
        raw_documents : iterable
            an iterable which yields either str, unicode or file objects
        Returns
        -------
        self : LogEntropyVectorizer
        """
        X = super(LogEntropyVectorizer, self).fit_transform(raw_documents)

        n_samples, n_features = X.shape
        gf = np.ravel(X.sum(axis=0)) # count total number of each words

        if self.smooth_idf:
            n_samples += int(self.smooth_idf)
            gf += int(self.smooth_idf)

        P = (X * sp.spdiags(1./gf, diags=0, m=n_features, n=n_features)) # probability of word occurence
        p = P.data
        P.data = (p * np.log2(p) / np.log2(n_samples))
        g = 1 + np.ravel(P.sum(axis=0))
        f = np.log2(1 + X.data)
        X.data = f
        # global weights
        self._G = sp.spdiags(g, diags=0, m=n_features, n=n_features)
        return self


    def fit_transform(self, raw_documents, y=None):
        self.fit(raw_documents)
        return self.transform(raw_documents)


    def transform(self, raw_documents):
        X = super(LogEntropyVectorizer, self).transform(raw_documents)
        check_is_fitted(self, '_G', 'global weight vector is not fitted')
        L = X * self._G  # sparse entropy matrix

        if self.norm is not None:
            L = normalize(L, norm=self.norm, copy=False)
        return L