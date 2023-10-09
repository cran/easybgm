# easybgm works for bgms

    Code
      summary(res_bgms)
    Output
      
       BAYESIAN ANALYSIS OF NETWORKS 
       Model type: ordinal 
       Number of nodes: 5 
       Fitting Package: bgms 
      --- 
       EDGE SPECIFIC OVERVIEW 
                Relation Estimate Posterior Incl. Prob. Inclusion BF     Category
        intrusion-dreams    0.472                 1.000          Inf     included
         intrusion-flash    0.587                 1.000          Inf     included
            dreams-flash    0.273                 0.888        7.929 inconclusive
         intrusion-upset    0.002                 0.088        0.096     excluded
            dreams-upset    0.156                 0.613        1.584 inconclusive
             flash-upset    0.102                 0.474        0.901 inconclusive
       intrusion-physior    0.136                 0.586        1.415 inconclusive
          dreams-physior    0.000                 0.045        0.047     excluded
           flash-physior    0.148                 0.606        1.538 inconclusive
           upset-physior    0.774                 1.000          Inf     included
      
       Bayes Factors larger than 10 were considered sufficient evidence for the classification 
       --- 
       AGGREGATED EDGE OVERVIEW 
       Number of edges with sufficient evidence for inclusion: 3 
       Number of edges with insufficient evidence: 5 
       Number of edges with sufficient evidence for exclusion: 2 
       Number of possible edges: 10 
       
       --- 
       STRUCTURE OVERVIEW 
       Number of visited structures: 54 
       Number of possible structures: 1024 
       Posterior probability of most likely structure: 0.156 
      ---

# easybgm works for bdgraph

    Code
      summary(res_bdgraph)
    Output
      
       BAYESIAN ANALYSIS OF NETWORKS 
       Model type: ggm 
       Number of nodes: 5 
       Fitting Package: bdgraph 
      --- 
       EDGE SPECIFIC OVERVIEW 
                Relation Estimate Posterior Incl. Prob. Inclusion BF     Category
        intrusion-dreams    0.379                  1.00          Inf     included
         intrusion-flash    0.439                  1.00          Inf     included
            dreams-flash    0.262                  1.00          Inf     included
         intrusion-upset    0.277                  1.00          Inf     included
            dreams-upset    0.182                  1.00          Inf     included
             flash-upset    0.020                  0.21        0.266 inconclusive
       intrusion-physior    0.000                  0.00        0.000     excluded
          dreams-physior    0.013                  0.27        0.370 inconclusive
           flash-physior    0.285                  1.00          Inf     included
           upset-physior    0.595                  1.00          Inf     included
      
       Bayes Factors larger than 10 were considered sufficient evidence for the classification 
       --- 
       AGGREGATED EDGE OVERVIEW 
       Number of edges with sufficient evidence for inclusion: 7 
       Number of edges with insufficient evidence: 2 
       Number of edges with sufficient evidence for exclusion: 1 
       Number of possible edges: 10 
       
       --- 
       STRUCTURE OVERVIEW 
       Number of visited structures: 4 
       Number of possible structures: 1024 
       Posterior probability of most likely structure: 0.592 
      ---

