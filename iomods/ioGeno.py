#!/usr/bin/env python

# Copyright Harald Grove, CIGENE, 2009
# This is GNU GPL Software: http://www.gnu.org/

# Description:
# A library with some genotype functions

import sys
import libPed
import libMark
import gzip

class Geno(object):
    def __init__(self, infile, pedfile=None, markfile=None):
        zipfile = False
        self.infile = infile
        self.pedfile = pedfile
        self.markfile = markfile
        if infile[-3:] == '.gz': zipfile = True
        self.genotypes = {} # Genotypes are stored as lists with animalID as key
        if zipfile: self.fgeno = gzip.open(infile,'r')
        else: self.fgeno = open(infile,'r')
        self.mafinfo = {}
        self.importFile()
        self.fgeno.close()

    def __getitem__(self, key):
        animal,marker = key[0],self.mark[str(key[1])][1] # Input is the markerID as found in the markerobject
        return self.genotypes[animal][marker*2:marker*2+2]

    def __setitem__(self,key,value):
        if value[0] not in ['0','1','2','3','4','5','A','C','G','T'] or value[1] not in ['0','1','2','3','4','5','A','C','G','T']:
            print "Setting illegal genotype",key,value
            sys.exit(1)
        try:
            animal,marker = key[0],self.mark[key[1]][1] # Input is the markerID as found in the markerobject
            self.genotypes[animal][marker*2:marker*2+2] = value
        except (KeyError,ValueError):
            print "Not recognized key",key
            sys.exit(1)

    def importFile(self):
        """ Merges duplicate samples """
        self.ped = libPed.Ped()
        firstline = True
        for line in self.fgeno:
            if line.strip().startswith('#'):
                if firstline:
                    line = line.strip('#').strip().split()
                    self.mark = libMark.Mark(line)
                    firstline = False
                continue
            l = line.strip().split()
            if len(l) == len(self.mark)+1 or len(l) == len(self.mark)*2+1: animal,sire,dam,geno = l[0],'0','0',l[1:]
            elif len(l) == len(self.mark)+3 or len(l) == len(self.mark)*2+3: animal,sire,dam,geno = l[0],l[1],l[2],l[3:]
            elif len(l) == 0: pass
            else:
                sys.stdout.write('Found %d genotypes for %s, expected %d markers\n' % (len(l)-3,l[0],len(self.mark)) )
                continue
                #sys.exit(1)
            newgeno = []
            # This part is to fix files where the 2 alleles have been combined into one
            if len(geno) == len(self.mark):
                for a in geno:
                    if len(a) == 2:
                        newgeno += [trans(a[0])] + [trans(a[1])]
                    elif len(a) == 1:
                        newgeno += [trans(a)] + [trans(a)]
                    elif len(a) == 4: # Includes either 'DEL' or 'INS'
                        if a[:3] == 'DEL': newgeno += ['D'] + [trans(a[3])]
                        elif a[1:] == 'DEL': newgeno += [trans(a[0])] + ['D']
                        elif a[:3] == 'INS': newgeno += ['I'] + [trans(a[3])]
                        elif a[1:] == 'INS': newgeno += [trans(a[0])] + ['I']
                        else: 
                            sys.stderr.write('Error in importing, unknown allele: %s\n' % a )
                            sys.exit(1)
                    elif len(a) == 6:
                        if a == 'DELDEL': newgeno += ['D'] + ['D']
                        elif a == 'DELINS': newgeno += ['D'] + ['I']
                        elif a == 'INSDEL': newgeno += ['I'] + ['D']
                        elif a == 'INSINS': newgeno += ['I'] + ['I']
                        else:
                            sys.stderr.write('Error in importing, unknown allele: %s\n' % a )
                            sys.exit(1)
                        
                geno = newgeno
            self.ped.addAnimal(animal,dam,sire,'F0','3')
            self.addGenotype(animal,geno)
            self.ped.updateSex()
        

    def addGenotype(self,animal,genos):
        """ Used for adding new animals to the data, genos is a list of genotypes
        """
        corr1 = corr2 = corr3 = 0
        if animal in self.genotypes: # used when two animals with the same name are present in the data
            for i in xrange(len(self.mark)):
                oa1,oa2 = self.genotypes[animal][2*i:2*i+2]
                na1,na2 = genos[2*i:2*i+2]
                if [na1,na2] == [oa1,oa2] or [na2,na1] == [oa1,oa2]: continue
                if oa1 == oa2 == '0':
                    self.genotypes[animal][2*i:2*i+2] = [na1,na2]
                    corr1 += 1
                elif na1 == na2 == '0': corr2 += 1
                else:
                    self.genotypes[animal][2*i:2*i+2] = ['0','0'] # If it's not possible to fill in later, then don't guess here
                    corr3 += 1
        else:
            self.genotypes[animal] = genos
        if (corr1+corr2+corr3) > 0: print "Non-identical genotypes in duplicated samples: %d in %s" % (corr3, animal)


    def updateMarkers(self):
        """ Updates the markerinfo with allele information
        """
        for mark in self.mark.getMarkers():
            a1 = a2 = '0'
            for animal in self.genotypes:
                if a1 != '0' and a2 != '0': break
                try: allele = self.__getitem__([animal,mark])
                except (KeyError,IndexError): continue
                if len(allele) != 2:
                    print "updateMarkers: not enough genotypes for %s and %s [%s]" % ( mark,animal,str(allele) )
                    sys.exit(1)
                if '0' in allele: continue
                try:
                    if a1 == '0': a1 = allele[0]
                except IndexError:
                    print "ioGeno, updateMarkers error: %s %s %s %d %d" % ( animal,mark,str(allele),len(self.genotypes[animal]), len(self.mark.getMarkers()) )
                    sys.exit(1)
                if a2 == '0' and a1 != '0': 
                    if a1 != allele[0]: a2 = allele[0]
                    elif a1 != allele[1]: a2 = allele[1]
            if a2 == '0': a2 = a1 # Set marker to homozygous if only one allele was found
            if a1 not in ['0','1','2','3','4','5','9'] or a2 not in ['0','1','2','3','4','5','9']:
                pass
                #print "<updateMarker>",a1,a2,mark,animal
            self.mark['alleles',mark] = sorted([a1,a2])
    
    def getMarkers(self):        
        """
        """
        self.updateMarkers()
        return self.mark
                
    def getPedigree(self):
        """
        """
        return self.ped
        
    def getPed(self):
        return self.ped
    
    def getAnimals(self):
        """
        """
        return [anim for anim in self.genotypes]
    
    def getGenotype(self,animal):
        """
        """
        return self.genotypes[animal]
        
    def getHaplotype(self,animal):
        """
        """
        return self.genotypes[animal][0::2][:],self.genotypes[animal][1::2][:]
        
    def getMAF(self,marker):
        try: return self.mafinfo[marker]
        except: return self.getMAF_(marker)
    
    def getMAF_(self,marker):
        """ Returns [minor allele,major allele, maf]
        """
        m1,m2,c1,tot = '0','0',0,0
        for animal in self.ped.getAnimals():
            a1,a2 = self.__getitem__([animal,marker])
            if '0' in a1+a2 or '-' in a1+a2: continue
            tot += 2
            if m1 == '0': m1 = a1
            if m2 == '0' and m1 != a1: m2 = a1
            elif m2 == '0' and m1 != a2: m2 = a2
            if a1 == m1: c1 += 1
            if a2 == m1: c1 += 1
        if tot == 0: return '0','0',0
        if c1 > tot/2.0:
            self.mafinfo[marker] = [m2,m1,1-c1/float(tot)]
            return m2,m1,1-c1/float(tot)
        else:
            self.mafinfo[marker] = [m1,m2,c1/float(tot)]
            return m1,m2,c1/float(tot)
            
    def setGenotype(self,animal,geno):
        if animal in self.genotypes: self.genotypes[animal] = geno
            
    def __str__(self):
        sep = '\t'
        out = '#'+sep+(sep+sep).join(self.mark.getMarkers())+'\n'
        for animal in self.ped.getAnimals():
            out += animal+sep+self.ped.getSire(animal)+sep+self.ped.getDam(animal)+sep+sep.join(self.genotypes[animal])+'\n'
        return out
        
    def fastPrint(self,fout):
        sep = '\t'
        fout.write('#'+sep+(sep+sep).join(self.mark.getMarkers())+'\n')
        for animal in self.ped.getAnimals():
            fout.write(animal+sep+self.ped.getSire(animal)+sep+self.ped.getDam(animal)+sep+sep.join(self.genotypes[animal])+'\n')
            
        
    def printData(self, fout, ped = [],markers = []):
        """ Use when modifications have been done to the markerlist
        """
        sep = '\t'
        if ped == []: ped = self.ped
        if markers == []: markers = self.mark
        fout.write('#'+sep+sep.join(markers.getMarkers())+'\n')
        for animal in ped.getAnimals():
            fout.write(animal+sep+ped.getSire(animal)+sep+ped.getDam(animal))
            fout.write(sep+sep.join([sep.join(self.__getitem__([animal,mark])) for mark in markers.getMarkers()])+'\n')

def writeData(outfile, cm, ped, markobj, chrom, warn = False, outfile2 = None, opt = '', pedfilter = None, markfilter=None):
    """ Writes out new official CIGENE genotype format
        warn = False means it will stop if trying to access markers or animals not present
    """
    sep = '\t'
    warned = False
    fout = open(outfile,'w')
    markers = markobj.getMarkers(chrom)
    fout.write('#'+sep+'\t\t'.join(markers)+'\n')
    for animal in ped.getAnimals():
        fout.write(animal+sep+ped.getSire(animal)+sep+ped.getDam(animal))
        for mark in markers:
            try: a1,a2 = cm[animal,mark]
            except KeyError:
                if warn:
                    sys.stderr.write('Unknown animal %s or marker %s, aborting...\n' % (animal,mark) )
                    sys.exit(1)
                else:
                    a1,a2 = '0','0'
                    if not warned: # only writes this warning once
                        sys.stdout.write('WARNING: Unknown animal(%s) or marker(%s)\n' % (animal,mark))
                        warned = True
            fout.write(sep+trans(a1)+sep+trans(a2))
        fout.write('\n')
    fout.close()

def trans(c):
    if c in 'A1': return '1'
    elif c in 'C2': return '2'
    elif c in 'G3': return '3'
    elif c in 'T4': return '4'
    elif c in 'D5': return '5'
    elif c in 'I6': return '6'
    else: return '0'
 
def main():
    print "Not a standalone program, exiting."
    
class EmptyGeno(object):
    def __init__(self):
        """ 
        """
        pass

    def __getitem__(self, key):
        return ['0','0']

    def getMAlleles(self, marker):
        return '-1','-1'
            
    def __str__(self):
        return ''
        
    def getAlleles(self,animal,marker,opt = 'lett'):
        return '0','0'
        
if __name__ == '__main__':
    #import cProfile
    #cProfile.run('main()')
    main()
