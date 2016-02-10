#! /usr/bin/env python
# Main input: Phased genotypes
#             Genotypes for new markers to place
# Reads pedigree and marker information from the phased genotypes file
#
# Harald Grove (c) 2010

import sys
import argparse
import numpy as np


class BuildMap(object):
    """
    """

    def __init__(self, options):
        self.ped = {}
        self.offspring = []
        self.readPedigree(options.phased)
        self.gen1, self.mark1 = self.readGenotypes(options.phased)
        self.gen2, self.mark2 = self.readGenotypes(options.newsnp)
        self.flog = open(options.logfile, 'w')
        self.outmat = options.outdist
        self.outlist = options.outmark
        self.phases = {0:np.zeros((len(ped),len(self.mark1)))}
        self.v = options.verbose
        
    def readGenotypes(self, infile):
        """
        Reads in genotypes from a Genotype file
        :param infile: Genotype file, first line must be '# mark1 mark2 ...'
        :return: gen1, genotypes as numpy; mark1, markers as list
        """
        fin = open(infile,'r')
        row = 0
        gen1 = {}
        for line in fin:
            if line.startswith('#'):
                l = line.strip().split()
                mark1 = l[1:]
                continue
            l = line.strip().split()
            gen1[l[0]] = [l[g]+l[l[g+1]] for g in range(3,len(l),2)]
            row += 1
        try:
            len(mark1)
        except:
            raise Exception("Missing header line in file %s" % infile)
        return gen1, mark1

    def readPedigree(self, infile):
        """
        Reads a genotype file and extract the pedigree from the first 3 columns
        :param infile:
        :return: None, populates self.ped
        """
        fin = open(infile,'r')
        rank = 0
        for line in fin:
            if line.startswith('#'):
                continue
            l = line.strip().split(None,3)
            try:
                name, father, mother = l[0], l[1], l[2]
            except IndexError:
                continue
            self.ped[name] = {'father':father,'mother':mother,'rank':rank}
            if father != '0' or mother != '0':
                self.offspring.append(name)
            rank += 1

    def findPhase(self, animal=''):
        """
        Schedules all or one animal for phase detection
        :param animal: the animal to generate phase for or leave empty for phasing all
        :return: self.phases, populates with phaseinfo
        """
        if animal == '':
            for animal in self.offspring:
                father, mother = self.ped[animal]['father'], self.ped[animal]['mother']
                self.setPhase(animal, father, mother)
        else:
            father, mother = self.ped[animal]['father'], self.ped[animal]['mother']
            self.setPhase(animal, father, mother)

    def setPhase(self, animal, father, mother):
        """

        :param animal:
        :param father:
        :param mother:
        :return:
        """
        row = self.ped[animal]['rank']
        angen = self.gen1[animal]
        fagen = self.gen1[father]
        mogen = self.gen1[mother]
        for i,ga in enumerate(angen):
            fa = fagen[i]
            ma = mogen[i]
            if fa[0] != fa[1]:
                if ga[0] == fa[0]:
                    self.phases[0][row,i] = -1
                elif ga[0] == fa[1]:
                    self.phases[0][row,i] = 1
                else:
                    raise Exception("Mendelian error between offspring %s and father %s [%s,%s] != [%s]"
                                    % (animal,father,ga[0],ga[1],fa[0]))
            else:
                self.phases[0][row,i] = 0
            if ma[0] != ma[1]:
                if ga[0] == ma[0]:
                    self.phases[1][row,i] = -1
                elif ga[0] == ma[1]:
                    self.phases[1][row,i] = 1
                else:
                    raise Exception("Mendelian error between offspring %s and mother %s [%s,%s] != [%s]"
                                    % (animal,mother,ga[0],ga[1],ma[1]))
            else:
                self.phases[1][row,i] = 0
            


    def impPhase(self, animal):
        """ Expands phaseinformation to include uninformative sites
        """
        if self.v: print
        "impPhase %s" % animal
        markers = self.oldlist
        self.impPhases[animal] = [[], []]
        sp, spos = self.findNextPhase(animal, 0, 0)
        dp, dpos = self.findNextPhase(animal, 1, 0)
        swind = [[sp, spos], [sp, spos]]  # previous,next
        dwind = [[dp, dpos], [dp, dpos]]  # previous,next
        for pos, mark in enumerate(markers):
            sph, dph = self.phases[animal][0][pos], self.phases[animal][1][pos]
            if sph != '-':
                self.impPhases[animal][0].append(sph)
                swind[0] = swind[1][:]
                sp, spos = self.findNextPhase(animal, 0, pos + 1)
                if spos == len(markers) and sp == '-':
                    swind[1] = [swind[0][0], len(markers)]
                else:
                    swind[1] = [sp, spos]
            elif sph == '-':
                if swind[0][0] == swind[1][0]:
                    self.impPhases[animal][0].append(swind[0][0])
                else:
                    self.impPhases[animal][0].append('-')
            if dph != '-':
                self.impPhases[animal][1].append(dph)
                dwind[0] = dwind[1][:]
                dp, dpos = self.findNextPhase(animal, 1, pos + 1)
                if dpos == len(markers) and dp == '-':
                    dwind[1] = [dwind[0][0], len(markers)]
                else:
                    dwind[1] = [dp, dpos]
            elif dph == '-':
                if dwind[0][0] == dwind[1][0]:
                    self.impPhases[animal][1].append(dwind[0][0])
                else:
                    self.impPhases[animal][1].append('-')

    def findNextPhase(self, animal, strand, start=0):
        """ Finds the phase in next informative position, including the indicated start-position """
        for pos in xrange(start, len(self.oldlist)):
            if self.phases[animal][strand][pos] != '-': return self.phases[animal][strand][pos], pos
        else:
            return '-', len(self.oldlist)

    def setHaplo(self):
        self.newgen = {}
        for nmark in self.newlist:
            for animal in self.ped.getAnimals():
                sire, dam = self.ped.getSire(animal), self.ped.getDam(animal)
                if sire == '0' and dam == '0': continue
                try:
                    if sire != '0': s1, s2 = self.c2[sire, nmark]
                except KeyError:
                    sire, s1, s2 = '0', '0', '0'
                try:
                    if dam != '0': d1, d2 = self.c2[dam, nmark]
                except KeyError:
                    dam, d1, d2 = '0', '0', '0'
                a1, a2 = self.c2[animal, nmark]
                if a1 == a2:
                    self.newgen[animal, nmark] = a1, a2
                elif sire != '0' and s1 == s2 and a1 == s1:
                    self.newgen[animal, nmark] = a1, a2
                elif sire != '0' and s1 == s2 and a2 == s1:
                    self.newgen[animal, nmark] = a2, a1
                elif dam != '0' and d1 == d2 and a1 == d1:
                    self.newgen[animal, nmark] = a2, a1
                elif dam != '0' and d1 == d2 and a2 == d1:
                    self.newgen[animal, nmark] = a1, a2
                else:
                    continue

    def setMarkerPhase(self):
        for pos, mark in enumerate(self.oldlist + ['dummy']):
            self.mphase[pos] = {}
            for anim in self.animlist:
                if self.ped.getSire(anim) == '0' and self.ped.getDam(anim) == '0':
                    # self.mphase[pos] += ['-','-']
                    continue
                if pos == 0:
                    self.mphase[pos][anim] = self.impPhases[anim][0][pos], self.impPhases[anim][1][pos]
                elif pos == len(self.oldlist):
                    self.mphase[pos][anim] = self.impPhases[anim][0][pos - 1], self.impPhases[anim][1][pos - 1]
                else:
                    p1, p2 = self.impPhases[anim][0][pos - 1], self.impPhases[anim][0][pos]
                    if p1 == p2:
                        ps = p1
                    elif p1 != p2:
                        ps = '-'
                    p1, p2 = self.impPhases[anim][1][pos - 1], self.impPhases[anim][1][pos]
                    if p1 == p2:
                        pd = p1
                    elif p1 != p2:
                        pd = '-'
                    self.mphase[pos][anim] = [ps, pd]

    def updateMarkerPhase(self, pos, mark):
        """ Updates the expected phase between each pair of markers
            Upgrade: Could be modified to only update what's needed
        """
        self.phases = {}
        self.findPhase()
        self.setMarkerPhase()

    # ********************************************************************************************

    def relocateMarkers(self, iterate=True):
        fmat = open(self.outmat, 'w')
        self.poslist = {}
        self.founders = self.ped.getFounders()
        for pos, nmark in enumerate(self.oldlist):
            fmat.write(nmark)
            self.poslist[nmark] = [999999, -1, 0]
            for gappos in xrange(0, len(self.mphase)):
                value = self.checkPlacement(nmark, gappos, pos)
                # print "%s, %d: %d" % (nmark,pos,value)
                if value[1] == 0:
                    fmat.write('\tNA')
                elif pos + 1 == gappos:
                    fmat.write('\t%.5f' % score)
                else:
                    score = (1.0 * value[0] / value[1])
                    fmat.write('\t%.5f' % score)
                    if score < self.poslist[nmark][0]:
                        self.poslist[nmark] = [score, gappos, 0]
                    elif score == self.poslist[nmark][0]:
                        if pos == gappos: self.poslist[nmark][1] = pos  # Prioritize old position if at minimum
                        self.poslist[nmark][2] = self.poslist[nmark][2] + 1
            fmat.write('\n')
            # Check if the marker has a new best position and if it should be moved
            score, newpos, shared = self.poslist[nmark]
            if iterate and score <= self.reclim and pos != newpos and nmark not in self.moved:  # Marker has a better position somewhere else
                self.moved[nmark] = 1
                self.flog.write('%s\t%d\t%d\n' % (nmark, pos, newpos))
                if newpos < pos:
                    self.oldlist.pop(pos)
                    self.oldlist.insert(newpos, nmark)
                elif newpos > pos:
                    self.oldlist.insert(newpos, nmark)
                    self.oldlist.pop(pos)
                fmat.close()
                return True
        fmat.close()
        return False

    def calcRecomb(self, options):
        flist = open(self.outlist, 'w')
        self.founders = self.ped.getFounders()
        for pos, nmark in enumerate(self.oldlist):
            flist.write(nmark)
            value = self.checkPlacement(nmark, pos, pos)
            score = (1.0 * value[0] / value[1])
            flist.write('\t%.5f\n' % score)
        flist.close()
        return False

    def checkPlacement(self, newmark, pos, oldpos):
        """ Counts number of recombinations by placing newmark in gap pos
            pos is the gap in front of the old marker in position pos
        """
        m1, m2 = self.newmark.getMAlleles(newmark)
        final = [0, 0]
        for parent in self.founders:
            p1, p2 = self.c2[parent, newmark]
            if p1 == p2: continue
            sex = self.ped.getSex(parent)
            strand = 1 - int(sex)
            count = [0, 0, 0, 0, 0]
            for off in self.ped.getOffspring(parent):
                current_p = self.mphase[pos][off][strand]
                if pos == oldpos:  # This is the oldposition
                    oldphase = self.impPhases[off][strand][oldpos]
                    if pos == 0:
                        p0 = self.mphase[pos + 1][off][strand]
                        if current_p == '-' and p0 == 'i':
                            current_p = 'o'
                        elif current_p == '-' and p0 == 'o':
                            current_p = 'i'
                    elif pos == len(self.oldlist):
                        p0 = self.mphase[pos + 1][off][strand]
                        if p0 == '-' and current_p == 'i':
                            current_p = 'o'
                        elif p0 == '-' and current_p == 'o':
                            current_p = 'i'
                    else:
                        p0 = self.mphase[pos + 1][off][strand]
                        if current_p == p0 == oldphase == '-':
                            continue
                        elif current_p != p0:
                            continue
                        elif current_p == p0 == '-':
                            if oldphase == 'o':
                                current_p = 'i'
                            elif oldphase == 'i':
                                current_p = 'o'
                if current_p == '-': continue
                try:
                    a = self.newgen[off, newmark]
                except KeyError:
                    continue
                a1 = a[strand]
                if a1 == '0': continue
                count[4] += 1
                if 0 < pos < len(self.oldlist):
                    recomb = 2
                else:
                    recomb = 1
                if current_p == 'i':
                    if a1 == m1:
                        count[0] += recomb
                    elif a1 == m2:
                        count[1] += recomb
                elif current_p == 'o':
                    if a1 == m1:
                        count[2] += recomb
                    elif a1 == m2:
                        count[3] += recomb
                else:
                    print
                "Missing if-statement in checkPosition"
            diff = min([count[0] + count[3], count[1] + count[2]])
            final[0] += diff
            final[1] += count[4]
        return final

    # ***************************************************************************************************

    def addMarkers(self, iterate=True):
        fmat = open(self.outmat, 'w')
        self.poslist = {}
        self.founders = self.ped.getFounders()
        for pos, nmark in enumerate(self.newlist):
            fmat.write(nmark)
            self.poslist[nmark] = [999999, -1, 0]
            for gappos in xrange(0, len(self.mphase)):
                value = self.checkPosition(nmark, gappos)
                # print "%s, %d: %d" % (nmark,pos,value)
                if value[1] == 0:
                    fmat.write('\tNA')
                else:
                    score = (1.0 * value[0] / value[1])
                    fmat.write('\t%.5f' % score)
                    if score < self.poslist[nmark][0]:
                        self.poslist[nmark] = [score, gappos, 0]
                    elif score == self.poslist[nmark][0]:
                        self.poslist[nmark][2] = self.poslist[nmark][2] + 1
            fmat.write('\n')
            # Check if the marker has a new best position and if it should be moved
            score, newpos, shared = self.poslist[nmark]
            if score <= self.reclim and iterate:  # Marker has a better position somewhere else
                self.flog.write('%s\t%d\n' % (nmark, newpos))
                self.oldlist.insert(newpos, nmark)
                self.updateMarkerPhase(newpos, nmark)
        fmat.close()
        return False

    def checkPosition(self, newmark, pos):
        """ Counts number of recombinations by placing newmark in gap pos
            pos is the gap in front of the marker in position pos
        """
        m1, m2 = self.newmark.getMAlleles(newmark)
        final = [0, 0]
        for parent in self.founders:
            p1, p2 = self.c2[parent, newmark]
            if p1 == p2: continue
            sex = self.ped.getSex(parent)
            strand = 1 - int(sex)
            count = [0, 0, 0, 0, 0]
            for off in self.ped.getOffspring(parent):
                current_p = self.mphase[pos][off][strand]
                if current_p == '-': continue
                try:
                    a = self.newgen[off, newmark]
                except KeyError:
                    continue
                a1 = a[strand]
                if a1 == '0': continue
                count[4] += 1
                if 0 < pos < len(self.oldlist):
                    recomb = 2
                else:
                    recomb = 1
                if current_p == 'i':
                    if a1 == m1:
                        count[0] += recomb
                    elif a1 == m2:
                        count[1] += recomb
                elif current_p == 'o':
                    if a1 == m1:
                        count[2] += recomb
                    elif a1 == m2:
                        count[3] += recomb
                else:
                    print
                "Missing if-statement in checkPosition"
            diff = min([count[0] + count[3], count[1] + count[2]])
            final[0] += diff
            final[1] += count[4]
        return final

    # *******************************************************************************************************
    # *******************************************************************************************************

    def createNewList(self, options):
        """ Writes a new markerfile
        """
        fout = open(self.outlist, 'w')
        if self.op == 'add':
            if not options.update: self.oldlist += ['dummy']
            for pos, mark in enumerate(self.oldlist):
                if not options.update:
                    for emark in self.poslist:
                        score, npos, dup = self.poslist[emark]
                        if score <= self.reclim and npos == pos: fout.write(
                            '%s\t%d\t%d\t%.5f\t%d\n' % (emark, pos, npos, score, dup))
                if mark in self.poslist:
                    score = self.poslist[mark][0]
                else:
                    score = 0
                if mark != 'dummy': fout.write('%s\t%d\tNA\t%.5f\tNA\n' % (mark, pos, score))
            for nmark in self.poslist:
                score, npos, dup = self.poslist[nmark]
                self.flog.write('%s\t%d\t%.5f\t%d\n' % (nmark, npos, score, dup))
        elif self.op == 're':
            for pos, mark in enumerate(self.oldlist):
                score, npos, dup = self.poslist[mark]
                fout.write('%s\t%d\t%d\t%.5f\t%d\n' % (mark, pos, npos, score, dup))
        fout.close()
        self.flog.close()

    def printPhase(self):
        """ Prints haplotypes to file
        """
        fout = open(self.out, 'w')
        if 0:
            for animal in self.ped.getAnimals():
                if animal not in self.phases: continue
                fout.write(animal + '\t' + 'P' + '\t' + ''.join(self.phases[animal][0]) + '\n')
                fout.write(animal + '\t' + 'M' + '\t' + ''.join(self.phases[animal][1]) + '\n')
        else:
            for animal in self.ped.getAnimals():
                if animal not in self.impPhases: continue
                fout.write(animal + '\t' + 'P' + '\t' + ''.join(self.impPhases[animal][0]) + '\n')
                fout.write(animal + '\t' + 'M' + '\t' + ''.join(self.impPhases[animal][1]) + '\n')
        fout.close()


def main():
    parser = argparse.ArgumentParser(description='Extends Linkage maps with new SNPs')
    parser.add_argument("-i", dest="phased", help="Phased input file", required=True)
    parser.add_argument("-j", dest="newsnp", help="Unplaced genotypes file")
    parser.add_argument("-o", dest="outmark", help="Final marklist after addition")
    parser.add_argument("-q", dest="outdist", help="Distance matrix")
    parser.add_argument("-l", dest="logfile", help="Log file")
    parser.add_argument("-e", dest="operation", help="Operation:[re/add/calc]", default="add")
    parser.add_argument("-v", dest="verbose", action="store_true", help="Prints runtime info", default=False)
    args = parser.parse_args()
    # Initiate the object, read in phased genotype data, new genotypes, pedigree
    imp = BuildMap(args)
    imp.findPhase()
    imp.setHaplo()
    if options.operation == 'add':
        imp.addMarkers(options.update)


if __name__ == '__main__':
    # import cProfile
    # cProfile.run('main()')
    main()
