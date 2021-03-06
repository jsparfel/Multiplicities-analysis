#!/usr/bin/python

# SUBMIT SEVERAL NICOLAS JOBS TO BATCH

# usage: submitJobs.py [-h] [-e EXEDIR] [-s EXESCRIPT] [-i INDIR] [-o OUTDIR]
#                      [-f] [-N] [-H] [-q]
#                      [periodFile]

# Submit several Nicolas jobs to batch

# positional arguments:
#   periodFile    D = "periods.dat"

# optional arguments:
#   -h, --help    show this help message and exit
#   -e EXEDIR     Exe dir name. D = "~/analyse/nicolas"
#   -s EXESCRIPT  Batch script "(analySIDIS|acceptance).csh". D =
#                 "analySIDIS.csh"
#   -i INDIR      Input dir name. D = "(trees|MCs)"
#   -o OUTDIR     Output dir name. D = "(rawmult|acceptance)"
#   -f, --force   Force even if output already exists
#   -N, --NOP     No OPeration
#   -H, --HELP    More help
#   -q, --quiet   Less verbosity

# NOTES:
# This is limited at present limited to:
# - analySIDIS_split and acceptance_split
# - 2016.
# Must be extended to other cases, viz.:
# - 2017, which would benefit a re-organisation of the input scheme, w/ e.g.
#  periods named <year><period> (15P09, 17P..)

# STATUS: "acceptance.csh" RUNNING ON ~ybedfer

from __future__ import print_function
import sys, os, argparse, re

exeScripts =   [ "analySIDIS.csh",   "acceptance.csh" ]
exeNames =     [ "analySIDIS_split", "accsplit" ]
inDirs =       [ "trees",            "MCs" ]
outDirs =      [ "rawmult",          "acceptance" ]
# Defaults
DPeriodFile = "periods.dat"
DExeDir = "~/analyse/nicolas"
DExeScript = exeScripts[0]
DInDir = inDirs[0]
DOutDir = outDirs[0]
# Lists of outputs
outList1 = ["DIS", "hadron", "hadron_theta", "xcheck", "DIS_zvtx", "hadron_pT", "hadron_zvtx", "count"]
outList2 = ["DIS/DIS", "DIS/DIS_zvtx", "hadron/hadron_pt", "hadron/hadron_theta", "hadron/hadron_zvtx", "electron/electron_pt", "electron/electron_theta", "electron/electron_zvtx"]
outLists = [ outList1, outList2 ]

##### PARSE COMMAND LINE
##### INPUT PERIODFILE
##### EXESCRIPT
##### INIT: INPUT-DIR, OUTPUT-DIR BASED ON EXESCRIPT
##### UPDATE INPUT-DIR, OUTPUT-DIR BASED ON ARGUMENTS
##### EXE-DIR: Check it contains exeName, exeScript, data
##### IN-DIR: Relative path ? prepend <EXEDIR>. Check it exists.
##### OUT-DIR: Relative path ? prepend <EXEDIR>. Check it exists.
#####          Update "2016" subdir and if required "2016/(DIS|hadron|electron)"
##### PARSE PERIOD FILE
#####   CHECK INPUT EXISTS
##### LOOP ON ACTIVE PERIODS AND SUBMIT
#####   CHECK OUTPUT DOES NOT ALREADY EXIST (unless force option)
#####   PASS ARG.S TO BATCH AS ENV. VAR.S

def printMoreHelp():
    print(' * submitJobs.py: Submit several Nicolas jobs to batch\n')

    print('- Periods to be processed specified via arg. [periodFile].')
    print('- One batch job per active period: <EXESCRIPT> =',end='')
    i = 0
    while i < len(exeScripts):
        print(' "%s"(%s)' % (exeScripts[i],exeNames[i]),end='')
        i += 1
    print()
    print('- Inputs TTree\'s specified in "<INDIR>/<period>/filelist.txt".')
    print('  Outputs go to "<OUTDIR>/2016/"...')
    print('   ...including "count.txt", renamed "counts_%d.txt,period", and the like.')
    print('  Relative <(IN|OUT)DIR> are relative to <EXEDIR> and not to $PWD.')

def main():

    ##### PARSE COMMAND LINE
    # Option MORE HELP
    if len(sys.argv) > 1 and (sys.argv[1] == '-H' or sys.argv[1] == '--HELP'):
        printMoreHelp()
        sys.exit(2)
    # "argparse" PARSER
    parser = argparse.ArgumentParser(description='Submit several Nicolas jobs to batch')
    parser.add_argument("periodFile",default=DPeriodFile,nargs='?',help="D = \"%s\"" % DPeriodFile)
    parser.add_argument('-e',dest='exeDir',default=DExeDir,help="Exe dir name. D = \"%s\"" % DExeDir)
    parser.add_argument('-s',dest='exeScript',default=DExeScript,help="Batch script \"(analySIDIS|acceptance).csh\". D = \"%s\"" % DExeScript)
    parser.add_argument('-i',dest='inDir',default=DInDir,help='Input  dir name. D = \"(%s|%s)\"' % (inDirs[0],inDirs[1]))
    parser.add_argument('-o',dest='outDir',default=DOutDir,help='Output dir name. D = \"(%s|%s)\"' % (outDirs[0],outDirs[1]))
    parser.add_argument('-f','--force',dest='force',action='store_true',help='Force even if output already exists')
    parser.add_argument('-N','--NOP',dest='nop',action='store_true',help='No OPeration')
    parser.add_argument('-H','--HELP',dest='moreHelp',action='store_true',help="More help")
    parser.add_argument('-q','--quiet',dest='quiet',action='store_true',help='Less verbosity')
    args = parser.parse_args()
    verbose = args.quiet == False
    force = args.force
    nop = args.nop

    ##### INPUT PERIODFILE
    # Open
    #periodFile = "local"
    # (use: "os.path.expanduser" to convert ~ into your home directory)
    periodFileName = os.path.expanduser(args.periodFile)
    try:
        periodFile = open(periodFileName,"r")
    except IOError as e:
        sys.stderr.write('** submitJobs: Error opening input <periodFile>: %s\n   => Exiting...\n' % str(e))
        sys.exit(2)

    ##### EXESCRIPT
    exeScript = args.exeScript
    ##### INIT: INPUT-DIR, OUTPUT-DIR BASED ON EXESCRIPT
    exeName = ""
    inDir = ""
    outDir = ""
    outputs = []
    i = 0
    match = False
    while i < len(exeScripts):
        if exeScript == exeScripts[i]:
            exeName = exeNames[i]
            inDir = inDirs[i]
            outDir = outDirs[i]
            outputs = outLists[i]
            match = True
            break
        i += 1
    if not match:
        sys.stderr.write('** submitJobs: Invalid <EXESCRIPT>(="%s") arg.\n   => Exiting...\n' % (exeScript))
        sys.exit(2)
    ##### UPDATE INPUT-DIR, OUTPUT-DIR BASED ON ARGUMENTS
    inDir = args.inDir
    outDir = args.outDir
    inconsistent = False
    if args.inDir == DInDir and exeScript != DExeScript:
        if exeScript == exeScripts[1]:
            inDir = inDirs[1]
        else:
            inconsistent = True
    if args.outDir == DOutDir and exeScript != DExeScript:
        if exeScript == exeScripts[1]:
            outDir = outDirs[1]
        else:
            inconsistent = True
    if inconsistent:
        sys.stderr.write('** submitJobs: Internal inconsistency.\n   => Exiting...\n')
        sys.exit(2)
        
    ##### EXE DIRECTORY: Check it contains exeScript, exeName, data
    m = re.match(r'^/',args.exeDir)
    if m:
        exeDir = args.exeDir
    else:
        exeDir = '%s/%s' % (os.environ['PWD'],args.exeDir)
    exeDir = os.path.expanduser(args.exeDir)
    exeScriptPath = "%s/%s" % (exeDir,exeScript)
    if not os.path.isfile(exeScriptPath):
        sys.stderr.write('** submitJobs: No "%s" file in <EXEDIR> = "%s"\n   => Exiting...\n' % (exeScript,exeDir))
        sys.exit(2)
    exePath = "%s/%s" % (exeDir,exeName)
    if not os.path.isfile(exePath):
        sys.stderr.write('** submitJobs: No "%s" file in <EXEDIR> = "%s"\n   => Exiting...\n' % (exeName,exeDir))
        sys.exit(2)
    dataPath = "%s/%s" % (exeDir,"data")
    if not os.path.isdir(dataPath):
        sys.stderr.write('** submitJobs: No "%s" sub-directory in <EXEDIR> = "%s"\n   => Exiting...\n' % ("data",exeDir))
        sys.exit(2)

    ##### IN-DIR: Relative path ? prepend <EXEDIR>. Check it exists.
    inDir = os.path.expanduser(inDir)
    m = re.match(r'^/',inDir)
    if m:
        inDirPath = inDir
        if not os.path.isdir(inDirPath):
            sys.stderr.write('** submitJobs: <INDIR> "%s" is not a directory\n   => Exiting...\n' % inDirPath)
            sys.exit(2)
    else:
        inDirPath = '%s/%s' % (exeDir,inDir)
        if not os.path.isdir(inDirPath):
            sys.stderr.write('** submitJobs: No <INDIR> "%s" sub-directory in <EXEDIR> = "%s"\n   => Exiting...\n' % (inDir,exeDir))
            sys.exit(2)
        
    ##### OUT-DIR: Relative path ? prepend <EXEDIR>. Check it exists.
    # Note: would be better to check output directory is writable
    outDir = os.path.expanduser(outDir)
    m = re.match(r'^/',outDir)
    if m:
        outDirPath = outDir
        if not os.path.isdir(outDirPath):
            sys.stderr.write('** submitJobs: <OUTDIR> "%s" is not a directory\n   => Exiting...\n' % outDirPath)
            sys.exit(2)
    else:
        outDirPath = "%s/%s" % (exeDir,outDir)
        if not os.path.isdir(outDirPath):
            sys.stderr.write('** submitJobs: No <OUTDIR> "%s" sub-directory in <EXEDIR> = "%s"\n   => Exiting...\n' % (outDir,exeDir))
            sys.exit(2)
        
    ##### OUT-DIR: Update "2016" subdir and if required "2016/(DIS|hadron|electron)"
    outDirPath = '%s/%s' % (outDirPath,"2016")
    if not os.path.isdir(outDirPath):
        print(' * submitJobs: Creating sub-directory \"%s\"' % outDirPath)
        if not nop:
            try:
                os.mkdir(outDirPath)
            except IOError as e:
                sys.stderr.write('** submitJobs: Error creating \"%s\": %s\n   => Exiting...\n' % (outDirPath,str(e)))
                sys.exit(2)
    if exeScript == "acceptance.csh":
        subDirs = [ "DIS", "hadron", "electron" ]
        for subDir in subDirs:
            outSubDir = "%s/%s" % (outDirPath,subDir)
            if not os.path.isdir(outSubDir):
                print(' * submitJobs: Creating sub-directory \"%s\"' % outSubDir)
                if not nop:
                    try:
                        os.mkdir(outSubDir)
                    except IOError as e:
                        sys.stderr.write('** submitJobs: Error creating \"%s\": %s\n   => Exiting...\n' % (outSubDir,str(e)))
                        sys.exit(2)   
        
    ##### PARSE PERIOD FILE
    if verbose:
        print(' * submitJobs: Parsing <perioFile> "%s" ...' % periodFileName)
    periods = []
    iLine = -1
    for line in periodFile:
        iLine += 1
        # Blank line: Skip
        m = re.match(r'\s*$',line)
        if m:
            continue
        # Comment: starts w/ # preceded by nothing else than \s (i.e.:
        # whitespace characters, which includes [ \t\n\r\f\v]).
        m = re.match(r'\s*#',line)
        if m:
            continue
        # Period entry: P[0-9][0-9] [01]
        m = re.match(r'\s*P\d\d\s+0',line)
        if m:  # Deactivated period
            continue
        m = re.match(r'\s*P\d\d\s+1',line)
        if m:  # Active period
            period = re.sub(r'\s*(P\d\d)\s+1\n$',r'\1',line)
            ##### CHECK INPUT EXISTS
            filelistDir = '%s/%s' % (inDirPath,period)
            if not os.path.isdir(filelistDir):
                sys.stderr.write('\n** submitJobs: Period "%s": No "%s" input dir\n   => Exiting...\n' % (period,filelistDir))
                sys.exit(2)
            filelistPath = '%s/filelist.txt' % filelistDir
            if not os.path.isfile(filelistPath):
                sys.stderr.write('\n** submitJobs: Period "%s": No "%s" input file\n   => Exiting...\n' % (period,filelistPath))
                sys.exit(2)
            periods.append(period)
            continue
        # Default: exit w/ error
        sys.stderr.write('\n** submitJobs: Unvalid entry in <periodFile> %s @ line %d: "%s"\n' % (periodFileName,iLine,line))
        sys.exit(2)
    periodFile.close()
    iLine += 1
    if not iLine:
        print('** submitJobs: No active period in <periodFile> "%s"\n   => Exiting...\n' % (periodFileName))
        sys.exit(2)
        
    ##### LOOP ON ACTIVE PERIODS AND SUBMIT
    print(' * submitJobs: Submitting \"%s\" jobs < \"%s\" > \"%s\"' % (exeScript,inDirPath,outDirPath))
    nActivePeriods = nSubmissions = 0
    for period in periods:
        nActivePeriods += 1
        exist = 0
        if not force:
            ##### CHECK OUTPUT DOES NOT ALREADY EXIST (unless force option)
            for output in outputs:  # Limiting ourselves to ".txt" files
                outFileName = '%s_%s.txt' % (output,period)
                outFilePath = '%s/%s' % (outDirPath,outFileName)
                if os.path.isfile(outFilePath):
                    exist = 1
                    print(' * submitJobs: Period "%s": Output file "%s" already exists in "%s"\n   => Skipping...\n' % (period,outFileName,outDirPath))
                    break
        if exist:
            continue
        print(' * submitJobs: Submitting period "%s"' % (period))
        nSubmissions += 1
        command = 'qsub -P P_compass'
        command += ' -V' # Submission shell env. passed to batch shell
        command += ' -l os=cl7'           # OS = CENTOS7
        # * Requested                                                          *
        # *   CPU cores:                    1 core(s)                          *
        # *   CPU time:                     02:46:40 (10000 seconds) (1)       *
        # * Consumed                                                           *
        # *   wallclock:                    00:15:56 (956 seconds)             *
        # *   CPU time:                     00:06:21 (381 seconds)             *
        # *   CPU scaling factor:           11.10                              *
        # *   normalized CPU time:          01:10:38 (4238 HS06 seconds)       *
        # *   CPU efficiency:               39 % (2)                           *
        command += ' -q long -l ct=10000' # CPU requested (see above)
        command += ' -l fsize=4G'         # Disk size requested
        # Book I/O on "/sps".
        # - Needed only if input TTree there.
        # - Could else be "-l xrootd=1" or "-l hpss=1"
        command += ' -l sps=1'
        command += ' -N aSs%s' % period   # Job name: a(naly)S(IDIS_)s(split)
        # Log file: <period>.log in outDir. Mixes stdout and stderr
        logFileName = '%s/%s.log' % (outDirPath,period)
        command += ' -o %s -j yes' % logFileName
        command += ' -V %s' % exeScriptPath
        if nop:   # NOP: Print submit command to stdout
            print('"%s"\n' % command)
            continue
        ##### PASS ARG.S TO BATCH AS ENV. VAR.S
        os.environ["analySIDIS_PERIOD"] = period
        os.environ["analySIDIS_EXEDIR"] = exeDir
        os.environ["analySIDIS_INPATH"] = inDirPath
        os.environ["analySIDIS_OUTPATH"] = outDirPath
        if verbose:
            print('"%s"\n' % command)
        status = os.system(command)
        if status:
            status = os.system(command)
            print(' * submitJobs: Try again submitting period "%s"' % period)
            if status:
                print('** submitJobs: Persistent submission error\n   => Exiting...\n')
                sys.exit(2)

    ##### SUMMARY
    print(' * submitJobs: %d periods submitted out of %d' % (nSubmissions,nActivePeriods))
    if nSubmissions == nActivePeriods:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
  main()
