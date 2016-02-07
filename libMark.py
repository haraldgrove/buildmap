#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2010
# This is GNU GPL Software: http://www.gnu.org/

# Description:
# A marker class
# Marker nr. is counted from 0

import sys
import random

class Mark(object):
    def __init__(self, infile= None):
        self.markers = {}
        self.alleles = {}
        self.marklist = []
        self.chromosomes = []
        if infile:
            try:
                self.fgeno = open(infile,'r')
                self.importFile()
                self.fgeno.close()
            except (IOError,TypeError):
                self.importFile(infile) # infile is a list of markers
        self.value = -1
        self.stop = len(self.marklist)-1
        self.offset = 0

    def importFile(self,ant = -1):
    
        def trans(s):
            if s in '1A': return '1'
            elif s in '2C': return '2'
            elif s in '3G': return '3'
            elif s in '4T': return '4'
            return s
        
        if ant == -1:
            countnr = 0
            style = ''
            warn = False
            for line in self.fgeno:
                if line.startswith('#'): continue
                l = line.strip().split()
                if len(l) == 0: continue
                # Figure out which marker-format the file is in
                if style == '':
                    if len(l) == 1: style = 'simple'
                    elif len(l) == 5 or (len(l) == 4 and len(l[2]) == 1 and len(l[3]) == 1): style = 'old'
                    else: style = 'new'
                if len(l) == 4 and style == 'new':
                    chrom,nr,name,pos = l
                    a1 = a2 = '0'
                elif style == 'old':
                    if len(l) == 4:
                        name,pos,a1,a2 = l
                        chrom = '99'
                    elif len(l) == 5: name,pos,a1,a2,chrom = l # attempts to use old marker file format
                    nr = countnr
                elif style == 'simple':
                    name = l[0]
                    pos = nr = countnr
                    chrom = '99'
                    a1 = a2 = '0'
                else:
                    sys.stderr.write('Incorrect marker format: [%s]\n' % (line.strip()))
                    sys.exit(1)
                self.markers[str(nr)] = [chrom,int(nr),name,pos]
                self.markers[name] = [chrom,int(nr),name,pos]
                self.marklist.append(name)
                self.chromosomes.append(chrom)
                if a1 != '0' and a2 != '0':
                    self.alleles[name] = trans(a1),trans(a2)
                    self.alleles[nr] = trans(a1),trans(a2)
                countnr += 1
            self.chromosomes = list(set(self.chromosomes))
        else:
            nr = 0
            for ma in ant:
                chrom,name,pos = '99',ma,nr
                self.markers[str(nr)] = [chrom,nr,name,pos]
                self.markers[name] = [chrom,nr,name,pos]
                self.marklist.append(name)
                nr += 1

    def addMarker(self,mark,chromosome,pos = -1):
        self.marklist.append(mark)
        if chromosome not in self.chromosomes: self.chromosomes.append(chromosome)
        nr = len(self.marklist)-1
        if pos == -1: pos = nr
        self.markers[mark] = [chromosome,int(nr),mark,pos]
        self.markers[str(nr)] = [chromosome,int(nr),mark,pos]
        
    def addAlleles(self,marker,a1,a2):
        self.alleles[marker] = trans(a1),trans(a2)
        self.alleles[str(self.markers[marker][1])] = trans(a1),trans(a2)

    def __getitem__(self, key):
        return self.markers[str(key)]
        
    def __setitem__(self,key,value):
        if 'alleles' in key:
            id,mark = key
            self.alleles[mark] = value
            self.alleles[self.markers[mark][1]] = value
        
    def getMarkers(self,chr='99'):
        if chr in ['99','']: return self.marklist
        return [m for m in self.marklist if self.markers[m][0] == chr]
        
    def getChroms(self):
        return self.chromosomes
        
    def getMAlleles(self,mark):
        return self.alleles[mark]
        
    def __iter__(self):
        return self
        
    def __len__(self):
        return len(self.marklist)
        
    def next(self):
        if self.value == self.stop:
            raise StopIteration
        self.value += 1
        return self.marklist[self.value]

    def __add__(self,other):
        """ Does not save any information beyond marker names and ranking """
        newM = self.marklist
        newM += [M2 for M2 in other.marklist if M2 not in self.markers]
        return Mark(newM)

    def __and__(self,other):
        newM = list(set(self.marklist).intersection(other.marklist))
        return Mark(newM)

    def __str__(self):
        out = ''
        sep = '\t'
        for mark in self.marklist:
            chr,nr,id,pos = self.markers[mark]
            out += chr+sep+str(nr)+sep+id+sep+str(pos)+'\n'
        return out
        
    def removeMarker(self,marker):
        pos = self.marklist.index(marker)
        self.marklist.pop(pos)

    def getOldMarkers(self):
        ut = []
        count = 1
        for marker in self.getMarkers():
            chr,nr,id,pos = self.__getitem__(marker)
            m1,m2 = self.getMAlleles(marker)
            ut.append([marker,str(pos),m1,m2,chr])
        return ut
        
    def writeSelection(self,r,outfile):
        """ Works only on bgl-style markerfiles
        """
        fout = open(outfile,'w')
        sep = '\t'
        for mark in self.marklist:
            if random.random()*100 > r: continue
            chr,nr,id,pos = self.markers[mark]
            try: m1,m2 = self.getMAlleles(mark)
            except:
                print "Requires bgl-style of markerfile"
                sys.exit(1)
            fout.write(id+sep+str(nr)+sep+m1+sep+m2+'\n')
        fout.close()

def main():
    try: m = Mark(sys.argv[1])
    except IOError:
        print "Error opening files"
        sys.exit(1)
    if sys.argv[2][:8] == 'rndsplit':
        m.writeSelection(int(sys.argv[2][8:]),sys.argv[3])
    else:
        print "Wrong option, use rndsplit#"

if __name__ == '__main__':
    main()
