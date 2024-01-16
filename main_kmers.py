#!/usr/bin/env python
# coding: utf-8

# In[3]:


class MainKmer( object ):
    """
    An object that holds the sequence of a main kmer, as well as a dictionary for each CG combo
    """
    def __init__(self, m_seq):
        self.seq = m_seq
        self.cg_combs = cg_split( self.seq )
        self.entropies = {}
    
    def calculate_freqs( self, sample ):
        """
        For a given sample, calculate the frequency of each possible cg combination 
        """
        counts = []
        total_count = 0
        
        for option in cg_combs:
            n = sample[option]
            counts.append( n )
            total_count += n
        
        # divide by the total number of counts to get the frequency
        freqs = []
        for c in counts:
            freqs.append( c / total_count )
        
        return freqs
    
    def calculate_entropy( self, sample, sample_name=None ):
        """
        Given a list of frequencies, calculate the Shannon entropy
        """
        import math
        
        # if no sample name is provided, create a sample name
        if not sample_name:
            sample_name = f"Sample{len(entropies)}"
        
        entropy = 0
        freqs = self.calculate_freqs( sample )
            
        # calculate the entropy
        for p in freqs:
            if( p != 0 ):
                entropy -= p * math.log(p, 2)
        
        # add the entropy to the dictionary
        if sample_name not in entropies:
            entropies[sample_name] = entropy
        
        return entropy
    


# In[4]:


def cg_split( seq ):
    """
    For a given sequence, return a list of all possible
    combinations of CG conversions
    """
    seq_len = len( seq )
    combos = []
    
    if( seq_len > 1 ):
        # if a CG is encountered, add a CG and a TG version
        if( seq[0:2] == 'CG'):
            for combo in cg_split( seq[2:] ):
                combos.append('CG' + combo)
                combos.append('TG' + combo)
        # otherwise, just move on to the next base
        else:
            for combo in cg_split( seq[1:] ):
                combos.append(seq[0] + combo)
    # if there are less than 2 bases left, just add the remaining bases
    else:
        combos = [ seq ]
    
    return combos


# In[23]:


if __name__ == "__main__":
    
    import pandas as pd
    
    # read in probe file and save to a dataframe
    probe_file = "Probes.csv"
    df = pd.read_csv(probe_file)
    
    # extract probe sequences
    probe_seqs = df['Sequence']
    
    
    kmer_size = 10 # can edit
    
    main_kmers = {} # empty dictionary to hold main kmer groups
    
    for probe in probe_seqs:
        seq = probe[8:] # trims off the /5Biosg/ prefix from sequence
        seq_len = len( seq ) # should be 120 bp
        
        count = 0
        for start_pos in range(0, seq_len - kmer_size + 1):
            count += 1
            # extract 10-mer sequence
            kmer = seq[start_pos : start_pos + kmer_size]
            
            # add it to the list of main kmers if not there already
            if kmer not in main_kmers:
                main_kmers[kmer] = MainKmer( kmer )
        
                    
    print(len(main_kmers)) # number of main kmers


# In[ ]:




