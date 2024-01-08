import sys
import os
import argparse
import random
import subprocess
from math import *

parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query",    type=str,        help="query", default="", required=False)
parser.add_argument("-s", "--status",   dest="status",   action="store_true")
parser.add_argument("-r", "--resubmit", dest="resubmit", action="store_true")
parser.add_argument("-k", "--kill",     dest="kill",     action="store_true")
args = parser.parse_args()

query = args.query
status = args.status
resubmit = args.resubmit
kill = args.kill

if query!= "":
   das_query = "dasgoclient --query='{QUERY}'".format(QUERY=query)
   print("DAS query: ",das_query)
   datasets = subprocess.check_output(das_query,shell=True)

   count = 0
   for ds in datasets.splitlines():
      ds = ((str(ds)).replace('b\'','')).replace('\'','')
      print(" ")
      print("--> Dataset: ",ds)
   
      ds_split = ds.split('/')
      request = str(ds_split[1]+"_"+ds_split[2])
      era = str(ds_split[2].split('_')[0])
   
      os.system('cp template_crab_cfg.py config_crab_'+str(count)+'.py') 
      with open('config_crab_'+str(count)+'.py') as f:
         newText=f.read().replace('REQUEST', str(request))
         newText=newText.replace('DATASET', str(ds)) 
         newText=newText.replace('ERA', str(era)) 
      with open('config_crab_'+str(count)+'.py', "w") as f:
         f.write(newText)
   
      submit_command = 'crab submit config_crab_'+str(count)+'.py'
      print("--> submit command: ",submit_command)
      os.system(submit_command)
      count += 1

elif status:
   dirs = subprocess.check_output('ls -d crab_*',shell=True)   
   for dir in dirs.splitlines():
      dir = ((str(dir)).replace('b\'','')).replace('\'','')
      print(" ")
      print("--> crab directory: ",dir)

      status_command = 'crab status '+str(dir)
      print("--> status command: ",status_command)
      os.system(status_command)

elif resubmit:
   dirs = subprocess.check_output('ls -d crab_*',shell=True)   
   for dir in dirs.splitlines():
      dir = ((str(dir)).replace('b\'','')).replace('\'','')
      print(" ")
      print("--> crab directory: ",dir)

      resubmit_command = 'crab resubmit '+str(dir)
      print("--> resubmit command: ",resubmit_command)
      os.system(resubmit_command)

elif kill:
   dirs = subprocess.check_output('ls -d crab_*',shell=True)   
   for dir in dirs.splitlines():
      dir = ((str(dir)).replace('b\'','')).replace('\'','')
      print(" ")
      print("--> crab directory: ",dir)

      kill_command = 'crab kill '+str(dir)
      print("--> kill command: ",kill_command)
      os.system(kill_command)
      os.system('rm -rf '+str(dir))


      

