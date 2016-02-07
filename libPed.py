#! /usr/bin/env python

# libPed v1.1
# Copyright 2009 Harald Grove

import sys

class Ped(object):
    
    def __init__(self,infile=None):
        self.families = {} # Based on animals
        self.animals = {} # fields: dam,sire,family,sex
        self.sires = {} # dict with offspring of sires
        self.dams = {} # dict with offspring of dams
        self.familyOrder = [] # List to record input order of families
        self.animalOrder = [] # Unique list of animals
        self.animalFullOrder = [] # Contains animals exactly as they are listed
        self.value = -1
        self.stop = -1
        try:
            self.fped = open(infile,'r')
            self.readPedigreeFile()
        except:
            pass
         
    def __getitem__(self,key):
        try: return self.animals[key]['ind']
        except: return key
            
    def __setitem__(self,key,value):
        id1,id2 = key
        self.animals[id2][id1] = value

    def getAnimals(self):
        """ Returns a list of animals in the order they were encountered """
        return self.animalOrder
        
    def getAllAnimals(self):
        """ Returns a list of every sample encountered in the input file """
        return self.animalFullOrder

    def getFamilies(self):
        """ Returns a list of families in the order they were encountered """
        return self.familyOrder
        
    def getFounders(self):
        return [anim for anim in self.animalOrder if self.animals[anim]['sire'] == self.animals[anim]['dam'] == '0']
        
    def getDam(self,animal):
        return self.animals[animal]['dam']
        
    def getSire(self,animal):
        return self.animals[animal]['sire']
    
    def getSex(self,animal):
        return self.animals[animal]['sex']
    
    def getFamily(self,animal):
        return self.animals[animal]['family']
    
    def getPhe(self,animal):
        return self.animals[animal]['phe']
    
    def getOffspring(self,value):
        if value in self.sires: return self.sires[value]
        elif value in self.dams: return self.dams[value]
        else: return []
        
    def getMembers(self,family):
        return self.families[family]
        
    def __len__(self):
        return len(self.animalOrder)
        
    def __iter__(self):
        return self
        
    def next(self):
        if self.value == self.stop: raise StopIteration
        self.value += 1
        return self.animalOrder[self.value]
    
    def __add__(self,other):
        newped = Ped()
        for anim in self.getAnimals():
            newped.addAnimal(anim,self.animals[anim]['dam'],self.animals[anim]['sire'],self.animals[anim]['family'][0],self.animals[anim]['sex'],self.animals[anim]['phe'])
        for anim in other.getAnimals():
            if anim not in self.animals:
                newped.addAnimal(anim,other.animals[anim]['dam'],other.animals[anim]['sire'],other.animals[anim]['family'][0],other.animals[anim]['sex'],self.animals[anim]['phe'])
        return newped
    
    def __radd__(self,other): 
        newped = Ped()
        for anim in self.getAnimals():
            newped.addAnimal(anim,self.animals[anim]['dam'],self.animals[anim]['sire'],self.animals[anim]['family'][0],self.animals[anim]['sex'],self.animals[anim]['phe'])
        for anim in other.getAnimals():
            if anim not in self.animals:
                newped.addAnimal(anim,other.animals[anim]['dam'],other.animals[anim]['sire'],other.animals[anim]['family'][0],other.animals[anim]['sex'],self.animals[anim]['phe'])
        return newped
        
    def __and__(self, other):
        newped = Ped()
        totanim = list(set(self.getAnimals()).intersection(other.getAnimals()))
        for anim in self.getAnimals():
            if anim in totanim:
                newped.addAnimal(anim,self.animals[anim]['dam'],self.animals[anim]['sire'],self.animals[anim]['family'][0],self.animals[anim]['sex'],self.animals[anim]['phe'])
        return newped
            
    def addAnimal(self,animal,dam,sire,family,sex='3',phen='-9'):
        """ A '0' for dam or sire indicates missing information
            If an animal is present more than once, merges the records provided they don't contradict
        """
        self.animalFullOrder.append([animal,sire,dam,family,sex,phen])
        if not animal in self.animals:
            self.animalOrder.append(animal)
            self.animals[animal] = {'dam':dam,'sire':sire,'family':[family],'sex':sex,'ind':len(self.animalOrder)-1,'phe':phen}
            if sire in self.animals: self.animals[sire]['sex'] = '1'
            if dam in self.animals: self.animals[dam]['sex'] = '0'
            self.stop += 1
        else:
            # Update the record of a previosly entered animal
            if self.animals[animal]['dam'] == '0': self.animals[animal]['dam'] = dam
            elif self.animals[animal]['dam'] == dam or dam == '0': pass
            else: 
                sys.stderr.write('%s is recorded with 2 different dams\n' % (animal))
                sys.exit(1)
            if self.animals[animal]['sire'] == '0': self.animals[animal]['sire'] = sire
            elif self.animals[animal]['sire'] == sire or sire == '0': pass
            else: 
                sys.stderr.write('%s is recorded with 2 different sires\n' % (animal))
                sys.exit(1)
            if family not in self.animals[animal]['family']:
                self.animals[animal]['family'].append(family)
            if self.animals[animal]['phe'] == '': self.animals[animal]['phe'] = phen
        try:
            if family in self.families: self.families[family].append(animal)
            elif family != '0':
                self.families[family] = [animal]
                self.familyOrder.append(family)
        except TypeError:
            print family
            print animal
            sys.exit(1)
        try: self.sires[sire].append(animal)
        except KeyError:
            if sire != '0': self.sires[sire] = [animal]
        try: self.dams[dam].append(animal)
        except KeyError:
            if dam != '0': self.dams[dam] = [animal]

    def updateSex(self):
        """ Assigns '1' and '0' to all known sires and dams, resp. """
        for sire in self.sires:
            if sire in self.animals: self.animals[sire]['sex'] = '1'
        for dam in self.dams:
            if dam in self.animals: self.animals[dam]['sex'] = '0'

    def anyDams(self):
        for animal in self.getAnimals():
            if self.getDam(animal) in self.animals: return True
        return False
        
    def readPedigreeFile(self):
        """ Reads pedigree as: animalID [sireID] [damID] [familyID] [sexID] """
        for line in self.fped:
            if line.startswith('#'): continue
            l = line.strip().split()
            animal = l[0]
            if len(l) > 1: sire = l[1]
            else: sire = '0'
            if len(l) > 2: dam = l[2]
            else: dam = '0'
            if len(l) > 3: family = l[3]
            else: family = 'F0'
            if len(l) > 4: sex = l[4]
            else: sex = '3'
            if len(l) > 5: phe = l[5]
            else: phe = ''
            self.addAnimal(animal, dam, sire,family,sex,phe)
        for sire in self.sires:
            if sire in self.animals: self.animals[sire]['sex'] = '1'
        for dam in self.dams:
            if dam in self.animals: self.animals[dam]['sex'] = '0'
    
    def printPedigree(self):
        sep = '\t'
        out = ''
        for family in self.getFamilies():
            for anim in self.getFamilyMembers(family):
                if anim in self.animals:
                    out += anim+sep+self.getSire(anim)+sep+self.getDam(anim)+sep+family+sep+self.getSex(anim)+'\n'
        return out

    def findGeneration(self):
        """ Assigns a generation number to each animal
        """
        founders = self.getFounders()
        self.generation = {}
        for parent in founders:
            gen = 1
            self.setGeneration(parent,gen)

    def setGeneration(self,parent,gen):
        if parent in self.generation:
            oldgen = self.generation[parent]
            if gen <= oldgen: return
        self.generation[parent] = gen
        offspring = self.getOffspring(parent)
        for off in offspring:
            self.setGeneration(off,gen+1)    

    def makeFamilies(self):
        """ Assigns unique family names to all related animals
        """
        self.families = {}
        self.familyOrder = []
        for animal in self.getAnimals(): self.animals[animal]['family'] = []
        familyID = 100
        self.assignedFamily = {}
        for animal in self.getAnimals():
            self.setFamilyName('F'+str(familyID),animal)
            familyID += 1

    def setFamilyName(self,famID,animal):
        if animal in self.assignedFamily or animal == '0': return
        self.assignedFamily[animal] = 1
        try: sire = self.getSire(animal)
        except KeyError: sire = '0'
        try: dam = self.getDam(animal)
        except KeyError: dam = '0'
        children = self.getOffspring(animal)
        self.setFamilyName(famID,sire)
        self.setFamilyName(famID,dam)
        for child in children:
            self.setFamilyName(famID,child)
        self.addAnimal(animal,'0','0',famID)

    def gather(self,animal,family,fout):
        """ Collects all animal that are connected in some way
        """
        if animal in self.printedAnimals or animal == '0': return # Already been here, or reached the end
        self.printedAnimals[animal] = 1
        sire,dam = self.getSire(animal),self.getDam(animal)
        children = self.getOffspring(animal)
        self.gather(sire,family,fout)
        self.gather(dam,family,fout)
        for child in children:
            self.gather(child,family,fout)
        fout.write(animal+'\t'+dam+'\t'+sire+'\t'+family+'\t'+self.getSex(animal)+'\n')
        
    def printFullped(self,outfile):
        """ Prints a pedigree gathering all connected animals into the same family
        """
        fout = open(outfile,'w')
        self.printedAnimals = {}
        families = self.getFamilies()
        for family in families:
            members = self.getFamilyMembers(family)
            out = self.gather(members[0],family,fout)
        fout.close()
        
    def __str__(self):
        sep = '\t'
        out = ''
        for family in self.families:
            members = self.getFamilyMembers(family)
            for animal in members:
                sire = self.getSire(animal)
                dam = self.getDam(animal)
                sex = self.getSex(animal)
                out += animal+sep+self.getSire(animal)+sep+self.getDam(animal)+sep+family+sep+self.getSex(animal)+'\n'
        return out
        
    def findKinship(self,outfile=None):
        """ Identifies all connections that are either
            identical
            siblings
            halfsibs
            sire - offspring
            grandsire - offspring
            uncle - offspring
            cousin - cousin
            Only works properly in bovine halfsibfamilies with unknown mother
        """
        try: fout = open(outfile,'w')
	except: pass
        kin = {}
        for animal1 in self.getAnimals():
            for animal2 in self.getAnimals():
                if animal1 == animal2:
                    kin[animal1,animal2] = 0 # Same
                    continue
                if (animal2,animal1) in kin: continue
                try: sire1 = self.getSire(animal1)
                except: sire1 = '0'
                try: sire2 = self.getSire(animal2)
                except: sire2 = '0'
                try: ff1 = self.getSire(sire1)
                except: ff1 = '0'
                try: ff2 = self.getSire(sire2)
                except: ff2 = '0'
                try: mf1 = self.getSire(self.getDam(animal1))
                except: mf1 = '0'
                try: mf2 = self.getSire(self.getDam(animal2))
                except: mf2 = '0'
                if (sire1 == animal2 and sire1 != '0') or (sire2 == animal1 and sire2 != '0'):
                    kin[animal1,animal2] = 1 # Sire
                elif (sire1 == sire2 and sire1 != '0'):
                    kin[animal1,animal2] = 2 # Half sib
                elif (ff1 == animal2 and ff1 != '0') or (ff2 == animal1 and ff2 != '0')  or (mf1 == animal2 and mf1 != '0') or (mf2 == animal1 and mf2 != '0'):
                    kin[animal1,animal2] = 3 # Grand sire
                elif (sire1 == ff2 and ff2 != '0') or (ff1 == sire2 and ff1 != '0') or (mf1 == sire2 and mf1 != '0') or (mf2 == sire1 and mf2 != '0'):
                    kin[animal1,animal2] = 4 # Uncle
                elif (ff1 == ff2 and ff1 != '0') or (ff1 == mf2 and ff1 != '0') or (ff2 == mf1 and ff2 != '0') or (mf1 == mf2 and mf1 != '0'):
                    kin[animal1,animal2] = 5 # Cousin
        if outfile:
            for pair in kin:
                animal1,animal2 = pair
                kinship = kin[pair]
                fout.write(animal1+'\t'+animal2+'\t'+str(kinship)+'\n')
            fout.close()
        return kin
    
    def findKinshipFull(self,outfile=None):
        """ Identifies all connections that are either
            identical
            siblings
            halfsibs
            sire - offspring
            grandsire - offspring
            uncle - offspring
            cousin - cousin
        """
        try: fout = open(outfile,'w')
        except: pass
        kin = {}
        for animal1 in self.getAnimals():
            for animal2 in self.getAnimals():
                if animal1 == animal2:
                    kin[animal1,animal2] = 0 # Same
                    continue
                if (animal2,animal1) in kin: continue
		#*********** Sire information ************
                try: s1 = self.getSire(animal1)
                except: s1 = None
                try: s2 = self.getSire(animal2)
                except: s2 = None
                try: ss1 = self.getSire(sire1)
                except: ss1 = None
                try: ss2 = self.getSire(sire2)
                except: ss2 = None
                try: sd1 = self.getSire(self.getDam(animal1))
                except: sd1 = None
                try: sd2 = self.getSire(self.getDam(animal2))
                except: sd2 = None
		#*********** Dam information *************
		try: d1 = self.getSire(animal1)
                except: d1 = None
                try: d2 = self.getSire(animal2)
                except: d2 = None
                try: ds1 = self.getSire(sire1)
                except: ds1 = None
                try: ds2 = self.getSire(sire2)
                except: ds2 = None
                try: dd1 = self.getSire(self.getDam(animal1))
                except: dd1 = None
                try: dd2 = self.getSire(self.getDam(animal2))
                except: dd2 = None
                if (s1 == animal2) or (s2 == animal1):
                    kin[animal1,animal2] = 1 # Sire
                elif s1 == s2 and d1 == d2 and s1 and d1:
                    kin[animal1,animal2] = 2 # Siblings
                elif (s1 == s2 and d1 != d2) or (d1 == d2 and s1 != s2) and s1 and d2:
                    kin[animal1,animal2] = 3 # Half-sibs
                #elif (ff1 == animal2 and ff1 != '0') or (ff2 == animal1 and ff2 != '0')  or (mf1 == animal2 and mf1 != '0') or (mf2 == animal1 and mf2 != '0'):
                #    kin[animal1,animal2] = 4 # Grand sire
                #elif (sire1 == ff2 and ff2 != '0') or (ff1 == sire2 and ff1 != '0') or (mf1 == sire2 and mf1 != '0') or (mf2 == sire1 and mf2 != '0'):
                #    kin[animal1,animal2] = 5 # Uncle
                #elif (ff1 == ff2 and ff1 != '0') or (ff1 == mf2 and ff1 != '0') or (ff2 == mf1 and ff2 != '0') or (mf1 == mf2 and mf1 != '0'):
                #    kin[animal1,animal2] = 6 # Cousin
                else:
                    kin[animal1,animal2] = 9
        if outfile:
            for pair in kin:
                animal1,animal2 = pair
                kinship = kin[pair]
                fout.write(animal1+'\t'+animal2+'\t'+str(kinship)+'\n')
            fout.close()
        return kin
   
    def clean(self):
        self.updateSex()
        animals = self.getAnimals()
        for animal in self.getAnimals():
            sire,dam = self.getSire(animal),self.getDam(animal)
            if sire not in animals: self.animals[animal]['sire'] = '0'
            if dam not in animals: self.animals[animal]['dam'] = '0'
        
 
    #***********************************************************
    # Old functions to keep compatibility
    #***********************************************************
    
    def getNumFamilies(self):
        return len(self.families)
        
    def getFamilyMembers(self,family):
        return self.families[family]
        
def main():
    ped = Ped(sys.argv[1])
    if sys.argv[2] == 'family':
        ped.makeFamilies()
        fout = open(sys.argv[3],'w')
        fout.write(str(ped))
        fout.close()
    elif sys.argv[2] == 'generation':
        ped.findGeneration()
        fout = open(sys.argv[3],'w')
        for anim in ped.generation:
            sire,dam,gen,fam,sex = ped.getSire(anim),ped.getDam(anim),ped.generation[anim],ped.getFamily(anim),ped.getSex(anim)
            famsize = len(ped.getMembers(fam[0]))
            fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (anim,sire,dam,fam[0],sex,gen,famsize) )
        fout.close()
    elif sys.argv[2] == 'clean':
        ped.clean()
        fout = open(sys.argv[3],'w')
        fout.write(str(ped))
        fout.close()
    else:
        print "Unknown option, use 'family' or 'generation'"

if __name__ == '__main__':
    main()
