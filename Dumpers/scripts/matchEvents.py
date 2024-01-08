#!/usr/bin/python

import argparse
import os
import math
import sys
from array import array
from ROOT import *

def getVal(tree,var):

  val = -999.
  if(var=='runId'): val = tree.runId
  if(var=='lumiId'): val = tree.lumiId
  if(var=='eventId'): val = tree.eventId     
  #else: print 'getVal ---> WARNING MISSING VAR: ',var
  #print("getVal: ",var,val)

  return val
  
if __name__ == '__main__': 

 gROOT.SetBatch(kTRUE)

 parser =  argparse.ArgumentParser(description='cloneTree')
 parser.add_argument('-t', '--inTree',    dest='inTree',    required=True,  type=str)
 parser.add_argument('-m', '--matchTree', dest='matchTree', required=True,  type=str) 
 args = parser.parse_args()
 
 inTree = args.inTree 
 matchTree = args.matchTree 

 treeMatch = TChain()
 treeMatch.AddFile(matchTree)
 print('Filling Ids from matchTree: ',treeMatch,' - nEntries=',treeMatch.GetEntries())
 
 matchIds = {}
 duplicateIds = {}
 for i,event in enumerate(treeMatch):

    if i>treeMatch.GetEntries():
      break 
    if i%100000==0:
      print("Reading Entry:",i)

    runId = getVal(treeMatch,'runId') 
    lumiId = getVal(treeMatch,'lumiId') 
    eventId = getVal(treeMatch,'eventId')
    id = tuple([int(runId),int(lumiId),int(eventId)])
       
    if id not in matchIds:
      matchIds[id] = 1
    else:
      #print("WARNING id=[",id,"] already in matchIds!")  
      matchIds[id] += 1
      
 
 treeIn = TChain()
 treeIn.AddFile(inTree) 
 print('Skimming inTrees: ',inTree,' - nEntries=',treeIn.GetEntries())

 outFile = TFile('Skimmed_tree.root','RECREATE')

 copyTree = treeIn.CloneTree(0);
 copyTree.SetName('bremdumper')
 copyTree.SetBranchStatus('*',1)
 for i,event in enumerate(treeIn):
   if i>treeIn.GetEntries():
     break 
   if i%100000==0:
     print("Reading Entry - ",i)   

   runId = getVal(treeIn,'runId') 
   lumiId = getVal(treeIn,'lumiId') 
   eventId = getVal(treeIn,'eventId')  
   id = tuple([int(runId),int(lumiId),int(eventId)])

   if id in matchIds: 
     #print("WARNING id=[",id,"] already in matchIds!")  
     continue

   if id not in duplicateIds:
     duplicateIds[id] = 1
   else:
     #print("WARNING id=[",id,"] already in duplicateIds!")  
     duplicateIds[id] += 1

   if duplicateIds[id]>1 :
     continue

   copyTree.Fill()

 outFile.Write()
 outFile.Close()

 
