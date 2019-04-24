class MyAlign:

    def __init__(self, lseqs, al_type = "protein"):
        self.listseqs = lseqs
        self.al_type = al_type
    
    def __len__(self): # number of columns
        return len(self.listseqs[0])
    
    def __getitem__(self, n):
        if type(n) is tuple and len(n) ==2: 
            i, j = n
            return self.listseqs[i][j]
        elif type(n) is int: return self.listseqs[n]
        return None
    
    def __str__(self):
        res = ""
        for seq in self.listseqs:
            res += "\n" + seq 
        return res
    
    def num_seqs(self):
        return len(self.listseqs)
   
       
    def column (self, indice):
        res = []
        for k in range(len(self.listseqs)):
            res.append(self.listseqs[k][indice])
        return res
    
    def consensus (self):
        # ...


if __name__ == "__main__": 
    alig = MyAlign(["ATGA-A","AA-AT-"], "dna")
    print(alig)
    print(len(alig))
    print(alig.column(2))
    print(alig[1,1])
    print(alig[0])
    print(alig.consensus())

    alig2 = MyAlign(["VJKK","JRSK","VRSK"])
    print(alig2.richpolarbasic())