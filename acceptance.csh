#!/bin/csh

# - Script to run acceptance on an arg. period
# Usage: qsub -V acceptance.csh
# - Arguments are passed as env. variables
#  (Note: names are taken from analySIDIS.)
#    - analySIDIS_PERIOD
#    - analySIDIS_EXEDIR
#    - analySIDIS_INPATH
#    - analySIDIS_OUTPATH
# - Works in combination w/ "./submitJobs.py".

##### CHECK ARG. ENV. VARIABLES analySIDIS_...
##### INSTALL LOCALLY I/O DIRS AND FILES AS REQUIRED BY acceptance_split
##### EXECUTE acceptance_split W/ OPTION -k -q
##### COPY OUTPUTS TO OUTPATH
#####   INCLUDING test.txt AND kinMC AND THE LIKE

printf '+++++++++++++++++++++++++++++++\n'
printf ' * acceptance: Run accsplit on an arg. period\n'
printf '+++++++++++++++++++++++++++++++\n'

##### WORKING DIR = $TMPDIR
if ( $?TMPDIR == 0 ) then
    printf '\n** acceptance: Inconsistency: No $TMPDIR\n'
else
    cd $TMPDIR
endif

##### CHECK ARG. ENV. VARIABLES analySIDIS_...
if ( $?analySIDIS_PERIOD == 0 ) then
    printf '\n** acceptance: Inconsistency: No $analySIDIS_PERIOD\n'
    exit 2
else
    set period = $analySIDIS_PERIOD
endif
set exeName = accsplit
if ( $?analySIDIS_EXEDIR == 0 ) then
    printf '\n** acceptance: Inconsistency: No $analySIDIS_EXEDIR\n'
    exit 2
else
    set exeDir = $analySIDIS_EXEDIR
    set exePath = $exeDir/$exeName
    if ( !( -e  $exePath ) ) then
	printf '\n** acceptance: Inconsistency: No "%s" in EXEDIR "%s"\n' $exeName $exeDir
	exit 2
    endif
endif
if ( $?analySIDIS_INPATH == 0 ) then
    printf '\n** acceptance: Inconsistency: No $analySIDIS_INPATH\n'
    exit 2
else
    set inPath = $analySIDIS_INPATH
endif
if ( $?analySIDIS_OUTPATH == 0 ) then
    printf '\n** acceptance: Inconsistency: No $analySIDIS_OUTPATH\n'
    exit 2
else
    set outPath = $analySIDIS_OUTPATH
endif

##### INSTALL LOCALLY I/O DIRS AND FILES AS REQUIRED BY accsplit
# PERIOD FILE
printf '%s 1\n' $period > periods.dat
printf '\n * Period file is:\n~~~~~~~~~~~\n'; cat periods.dat; printf '~~~~~~~~~~~\n'
# INPUT DIR
mkdir MCs
mkdir MCs/$period
\cp $inPath/$period/filelist.txt MCs/$period
# data DIR
# See "Inputs" in "acceptance_split.cc"
# Note:
#   We opt here for copying a selection of the contents of the "data"
#  sub-directory. Obviously it spares us from copying much unuseful stuff.
#   But is slightly dangerous, in that it's not foreward-compatible: i.e.
#  any future change to the "analySIDIS_split.cc" code implying new "data"
#  files would have to be ported to this script.
#   Yet as a safety precuation, in all such reading of inputs from
#  external files, the availibility of the files should be checked, and an
#  exception be raised whenever it's not granted. That's a deficiency of
#  the original source code of Nicolas that not checks were programmed and
#  the non availability of input files was passing unnoticed. But it has
#  been corrected since. If we stick to this safety rule, any future
#  addition to the set of inputs not propagated to the batch script, will
#  be detected.
#define target_file_2012 "data/target-107924-109081.dat"
#define target_file_2016 "data/target-274508-274901.dat"
mkdir data
set dataFiles = ( \
    target-107924-109081.dat target-274508-274901.dat )
foreach dataFile ( $dataFiles )
    cp $exeDir/data/$dataFile data
    if ( $status != 0 ) then
	printf '** acceptance: Error copying \"%s\" to $PWD\n' $exeDir/data/$dataFile
	exit 2
    endif
end

# OUTPUT DIR
mkdir acceptance
mkdir acceptance/2016
set dirs = ( DIS hadron electron )
foreach dir ( $dirs )
    mkdir acceptance/2016/$dir
end   

##### EXECUTE accsplit W/ OPTION -k -q
printf '\n~~~~~~~~~~~~~~~~~~~~~\n'
printf ' * Execute "%s"\n' $exePath
printf   '~~~~~~~~~~~~~~~~~~~~~\n'
$exePath periods.dat -k -q
set exeStatus = $status
printf   '~~~~~~~~~~~~~~~~~~~~~\n'
if ( $exeStatus != 0 ) then
    printf ' * acceptance: Execution returns %d\n' $exeStatus
printf   '~~~~~~~~~~~~~~~~~~~~~\n'
endif
    
##### COPY OUTPUTS TO OUTPATH
set copyStatus = 0
foreach dir ( $dirs )
    set files = `\ls acceptance/2016/$dir`
    printf ' * Copying to "%s":' $outPath/$dir; echo $files
    foreach file ( $files )
	\cp acceptance/2016/$dir/$file $outPath/$dir
	if ( $status != 0 ) then
	    set copyStatus = 1
	    printf '** Error copying "%s" to "%s"\n' $file $outPath/$dir
	endif
    end
end
#####   INCLUDING test.txt, kinMC.(pdf|root), Trigger_Coverage.pdf, trigmask.dat
set files = ( test.txt kinMC.pdf kinMC.root Trigger_Coverage.pdf trigmask.dat )
foreach file ( $files )
    if ( -e $file ) then
	if      ( $file == test.txt ) then
	    set fich = $outPath/`printf 'test_%s.txt' $period`
	else if ( $file == kinMC.pdf ) then
	    set fich = $outPath/`printf 'kinMC_%s.pdf' $period`
	else if ( $file == kinMC.root ) then
	    set fich = $outPath/`printf 'kinMC_%s.root' $period`
	else if ( $file == Trigger_Coverage.pdf ) then
	    set fich = $outPath/`printf 'Trigger_Coverage_%s.pdf' $period`
	else
	    set fich = $outPath/`printf 'trigmask_%s.dat' $period`
	endif
	printf ' * Copying "%s" to "%s"\n' $file $fich
	\cp $file $fich
	if ( $status != 0 ) then
	    set copyStatus = 1
	    printf '** Error copying "%s" to "%s"\n' $file $fich
	endif
    else
	if ( $exeStatus == 0 ) set exeStatus = 1
	printf '** Inconsistency: No output "%s"\n'
    endif
end

##### ERROR
set error = `expr $exeStatus \* 2 + $copyStatus`
printf '+++++++++++++++++++++++++++++++\n'
printf '\ * acceptance: Script execution status = %d\n' $error
printf '+++++++++++++++++++++++++++++++\n'
exit $copyStatus
