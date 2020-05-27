#! /usr/bin/python3

# Copyright Kyle Upton 2020. Not for any use or replication without express premission until i get a standard copyright applied.


import argparse, random, os, time
import datetime
import h5py

import numpy as np

from collections import OrderedDict as od
from timeit import default_timer as timer


def rev_comp(dna):
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV-', 'tgcayrkmvhdbTGCAYRKMVHDB-')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

class Chroms_setup:
    def __init__(self, faDir, chrs, gapFile):
        self.faDir = faDir
        self.chrs = chrs
        self.gapFile = gapFile
        self.chrSeq = od()
        self.gaps = {}
        self.lenDict = od()
        self.effLenDict = od()

    def sort_chroms(self, inDir, alts=False, rands=False):
        ''' sort chromosomes by intiger then string value'''
        files = [f for f in os.listdir(inDir) if f.endswith('.fa')]
        if not alts:
            files = [f for f in files if not 'alt' in f]
        if not rands:
            files = [f for f in files if not 'random' in f]
        strs = []
        ints = []
        for file in files:
            ID = file.lstrip('chr').rstrip('.fa')
            try:
                ints.append(int(ID))
            except ValueError:
                strs.append(ID)
        ints.sort()
        strs.sort()
        ints = [str(i) for i in ints]
        order = dict(zip((ints+strs), range(len(ints+strs))))
        files.sort(key = lambda x: order[x.lstrip('chr').rstrip('.fa')])
        return files
        
    def read_chrs(self):
        def process(file):
            with open(file, 'r') as f:
                lines = f.readlines()
                chrLen = 0
                chrCount = 0
                thisSeq = ''
                for line in lines:
                    line = line.strip()
                    if line.startswith('>'):
                        if chrCount > 0:
#                             if not thisChr.endswith('alt'):
                            self.lenDict[thisChr] = chrLen
                            self.chrSeq[thisChr] = thisSeq
                        thisChr = line[1:]
                        chrLen = 0
                        thisSeq = ''
                        chrCount += 1
                    else:
                        chrLen += len(line)
                        thisSeq += line
#                 if not thisChr.endswith('alt'):
                self.lenDict[thisChr] = chrLen
                self.chrSeq[thisChr] = thisSeq
                    
        if self.chrs==False:
            files = self.sort_chroms(self.faDir)
            print(files)
            for file in files:
                print(file)
                process(os.path.join(self.faDir,file))
        else:
            process(self.faDir + self.chrs)
        
    def find_gaps(self):
        gapDict = {}
        for chrom, seq in self.chrSeq.items():
            theseGaps = []
            c = 0
            gapStart = -1
            gapEnd = 0
            for base in seq:
                if ((base=='N') or (base=='n')):
                    if gapStart == -1:
                        gapStart = c
                        gapEnd = -1
                else:
                    if gapEnd == -1:
                        gapEnd = c
                        theseGaps.append([gapStart,gapEnd])
                        gapStart = -1
                        gapEnd = 0
                c+=1
            if gapEnd == -1:
                gapEnd = c
                theseGaps.append([gapStart,gapEnd])
            gapDict[chrom] = theseGaps
            
        gapFile = os.path.join(self.faDir, 'gaps.bed')
        with open(gapFile, 'w') as o:
            for chrom, gaps in gapDict.items():
                for gap in gaps:
                    o.write('\t'.join([chrom, str(gap[0]), str(gap[1])]) + '\n')
        return gapFile

    def read_gaps(self):
        if self.gapFile==False:
            self.gapFile = self.find_gaps()
        
        with open(self.gapFile, 'r') as bed:
            for line in bed:
                line = line.strip()
                if line.startswith('chr'):
                    # print line
                    data = line.split('\t')
                    # if not data[0] == 'chrY':
                    # if not data[0] == 'chrM':
                    if not data[0].endswith('alt'):
                        try:
                            self.gaps[data[0]].append([int(data[1]),int(data[2])])
                        except KeyError:
                            self.gaps[data[0]] = [[int(data[1]),int(data[2])]]
        return self.gapFile

    def calc_effective_sizes(self):
        self.effGapPosDict = od()
        self.effStartDict = {}
        self.effStartDict2 = {}
        self.effEndDict = {}
        effectiveStart = 0
        for chrom in self.lenDict.keys():
            self.effStartDict[chrom] = effectiveStart
            self.effStartDict2[effectiveStart] = chrom
            chromGapPos = []
            totalGaps = 0
            try:
                gaps = self.gaps[chrom]
            except KeyError:
                # print 'No entry for ' + chrom + ' in gaps dictionary'
                gaps = [[0,10]]
            for gap in gaps:
                chromGapPos.append(gap[0] - totalGaps)
                thisGap = gap[1] - (gap[0] - 1)
                totalGaps += thisGap
            fullLength = self.lenDict[chrom]
            effectiveLength = fullLength - totalGaps
            effectiveStart += effectiveLength
            self.effEndDict[chrom] = effectiveStart
            self.effLenDict[chrom] = effectiveLength
            self.effGapPosDict[chrom] = chromGapPos
        print (self.effLenDict)
        print (sum(self.effLenDict.values()))
        print (self.effStartDict)

    def chr_order(self):
        pass
        return False
        # Want a function to shuffle the chromosome order

    def setup_h5(self, start, h5Path, genome):
        if not os.path.exists(h5Path):
            os.makedirs(h5Path)
        with h5py.File(os.path.join(h5Path, genome + '.hdf5'), 'w') as f:
            self.chromSeqs      = f.create_group(genome + 'Chrs')
            self.chromGaps      = f.create_group(genome + 'Gaps')
            self.chromLens      = f.create_group(genome + 'Lens')
            self.chromEffLens   = f.create_group(genome + 'EffLens')
            self.chromEffStart  = f.create_group(genome + 'EffStart')
            self.chromEffEnd    = f.create_group(genome + 'EffEnd')
            self.chromEffGaps   = f.create_group(genome + 'EffGaps')
            for k,v in self.chrSeq.items():
                print ('starting np and h5 routine for chrom ' + k)
                start2 = timer()
                npSeq = np.array(list(v), dtype='|S1')
                stop = timer()
                print ('Populating npSeq took \t\t' + str(stop-start2) + ' seconds')
                self.chromSeqs.create_dataset(k, data=npSeq)
                stop = timer()
                print ('Seting up h5 dataset took \t' + str(stop-start2) + ' seconds')
                try:
                    npGap = np.array(list(self.gaps[k]), dtype=np.uint32)
                    print (self.gaps[k][:2])
                    print (self.gaps[k][-2:])
                except KeyError:
                    print ('Gap Key Error\t' + k)
                    npGap = np.array([0], dtype=np.uint32)
                self.chromGaps.create_dataset(k, data=npGap)
                for thisDict, thisDataSet in [
#                     [self.chrSeq,self.chromSeqs],
#                     [self.gaps,self.chromGaps],
                    [self.lenDict,self.chromLens],
                    [self.effLenDict,self.chromEffLens],
                    [self.effStartDict,self.chromEffStart],
                    [self.effEndDict,self.chromEffEnd],
                    [self.effGapPosDict,self.chromEffGaps]]:
                    try:
                        thisVar = np.array(thisDict[k], dtype=np.uint32)
                    except KeyError:
                        thisVar = np.array([0], dtype=np.uint32)
                    thisDataSet.create_dataset(k, data=thisVar)


def main(args):
	gapFile=os.path.join(args.inDir, 'gaps.bed')
	chromsSetup = Chroms_setup (args.inDir, args.chrs, gapFile)

	start = timer()
	chromsSetup.read_chrs()
	stop = timer()
	print ('reading chroms took ' + str(stop-start) + ' seconds')

	start = timer()
	chromsSetup.read_gaps()
	stop = timer()
	print ('reading gaps took ' + str(stop-start) + ' seconds')

	start = timer()
	chromsSetup.calc_effective_sizes()
	stop = timer()
	print ('calculating effective sizes took ' + str(stop-start) + ' seconds')

	start = timer()
	chromsSetup.setup_h5(start, args.outDir, args.genome)
	stop = timer()
	print ('Seting up h5 took ' + str(stop-start) + ' seconds')
	# chroms.chr_order()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load a genome into hdf5 format')

    parser.add_argument('-c', '--chrs', dest='chrs', default=False, help='absolute path for infut fasta file containing all chromosomes in sorted order')
    parser.add_argument('-i', '--indir', dest='inDir', default='~/fa/hg38', help='absolute path for infut fasta files')
    parser.add_argument('-o', '--outdir', dest='outDir', default='~/h5/', help='absolute path for output files in hdf5 format')
    parser.add_argument('-g', '--genome', dest='genome', default='hg38', help='prefix/name to be used for genome')

    args = parser.parse_args()
    main(args)
