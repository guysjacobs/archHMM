###Reads in a list of individuals separated by EOLs

###GSJ 05/10/2017

def readin_populationFileEOL(infile):
    with open(infile, 'rb') as f:
        pop = []
        for line in f:
            ind = line[0:-2] if line[-2:] == '\r\n' else line[0:-1] if (line[-1:] == '\r' or line[-1:] == '\n') else line
            pop.append(ind)
    return pop
