#!/bin/csh

# - Script to run analySIDIS on an arg. period
# Usage: qsub -V analySIDIS.csh
# - Arguments are passed as env. variables
#    - analySIDIS_PERIOD
#    - analySIDIS_EXEDIR
#    - analySIDIS_INPATH
#    - analySIDIS_OUTPATH
# - Works in combination w/ "./submitJobs.py".

##### CHECK ARG. ENV. VARIABLES analySIDIS_...
##### INSTALL LOCALLY I/O DIRS AND FILES AS REQUIRED BY analySIDIS_split
##### EXECUTE analySIDIS_split W/ OPTION -k -q
##### COPY OUTPUTS TO OUTPATH
#####   INCLUDING count.txt AND kinSIDIS

printf '+++++++++++++++++++++++++++++++\n'
printf ' * analySIDIS: Run analySIDIS_split on an arg. period\n'
printf '+++++++++++++++++++++++++++++++\n'

##### WORKING DIR = $TMPDIR
if ( $?TMPDIR == 0 ) then
    printf '\n** analySIDIS: Inconsistency: No $TMPDIR\n'
else
    cd $TMPDIR
endif

##### CHECK ARG. ENV. VARIABLES analySIDIS_...
if ( $?analySIDIS_PERIOD == 0 ) then
    printf '\n** analySIDIS: Inconsistency: No $analySIDIS_PERIOD\n'
    exit 2
else
    set period = $analySIDIS_PERIOD
endif
set exeName = analySIDIS_split
if ( $?analySIDIS_EXEDIR == 0 ) then
    printf '\n** analySIDIS: Inconsistency: No $analySIDIS_EXEDIR\n'
    exit 2
else
    set exeDir = $analySIDIS_EXEDIR
    set exePath = $exeDir/$exeName
    if ( !( -e  $exePath ) ) then
	printf '\n** analySIDIS: Inconsistency: No "%s" in EXEDIR "%s"\n' $exeName $exeDir
	exit 2
    endif
endif
if ( $?analySIDIS_INPATH == 0 ) then
    printf '\n** analySIDIS: Inconsistency: No $analySIDIS_INPATH\n'
    exit 2
else
    set inPath = $analySIDIS_INPATH
endif
if ( $?analySIDIS_OUTPATH == 0 ) then
    printf '\n** analySIDIS: Inconsistency: No $analySIDIS_OUTPATH\n'
    exit 2
else
    set outPath = $analySIDIS_OUTPATH
endif

##### INSTALL LOCALLY I/O DIRS AND FILES AS REQUIRED BY analySIDIS_split
# PERIOD FILE
printf '%s 1\n' $period > periods.dat
printf '\n * Period file is:\n~~~~~~~~~~~\n'; cat periods.dat; printf '~~~~~~~~~~~\n'
# INPUT DIR
mkdir trees
mkdir trees/$period
\cp $inPath/$period/filelist.txt trees/$period
# data DIR
# See "Inputs" in "analySIDIS_split.cc"
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
#define mat_RICH_2006_name "data/rich_mat_2006.txt"
#define mat_RICH_2016_name "data/rich_mat_2016.txt"
#define err_RICH_name "data/rich_mat_error.txt"
#define target_file_2012 "data/target-107924-109081.dat"
#define target_file_2016 "data/target-274508-274901.dat"
#define proton_sirc "data/proton_semi_inclusive_RC.txt"
#define proton_irc "data/hh160_r1998_f2tulay_compass_grv.asy_hcorr.txt"
#define ElectronPi "data/electron_pion_contamination.txt"
#define ElectronPiVtx "data/electron_pion_contamination_vtx.txt"
#define ElectronPiTheta "data/electron_pion_contamination_Theta.txt"
#define ElectronPipT "data/electron_pion_contamination_pT.txt"
mkdir data
set dataFiles = ( \
    rich_mat_2006.txt rich_mat_2016.txt rich_mat_2016_julien.txt rich_mat_error.txt rich_mat_error_julien.txt \
    target-107924-109081.dat target-274508-274901.dat \
    proton_semi_inclusive_RC.txt hh160_r1998_f2tulay_compass_grv.asy_hcorr.txt \
    electron_pion_contamination.txt electron_pion_contamination_vtx.txt \
    electron_pion_contamination_Theta.txt electron_pion_contamination_pT.txt )
foreach dataFile ( $dataFiles )
    cp $exeDir/data/$dataFile data
    if ( $status != 0 ) then
	printf '** analySIDIS: Error copying \"%s\" to $PWD\n' $exeDir/data/$dataFile
	exit 2
    endif
end

# OUTPUT DIR
mkdir rawmult
mkdir rawmult/2016/

##### EXECUTE analySIDIS_split W/ OPTION -k -q
printf '\n~~~~~~~~~~~~~~~~~~~~~\n'
printf ' * Execute "%s"\n' $exePath
printf   '~~~~~~~~~~~~~~~~~~~~~\n'
$exePath periods.dat -k -q
set exeStatus = $status
printf   '~~~~~~~~~~~~~~~~~~~~~\n'
if ( $exeStatus != 0 ) then
    printf ' * analySIDIS: Execution returns %d\n' $exeStatus
printf   '~~~~~~~~~~~~~~~~~~~~~\n'
endif
    
##### COPY OUTPUTS TO OUTPATH
set copyStatus = 0
set files = `\ls rawmult/2016/`
printf ' * Copying to "%s":' $outPath; echo $files
foreach file ( $files )
    \cp rawmult/2016/$file $outPath
    if ( $status != 0 ) then
	set copyStatus = 1
    	printf '** Error copying "%s" to "%s"\n' $file $outPath
    endif
end
#####   INCLUDING count.txt AND kinSIDIS.(pdf|root)
set files = ( count.txt kinSIDIS.pdf kinSIDIS.root )
foreach file ( $files )
    if ( -e $file ) then
	if      ( $file == count.txt ) then
	    set fich = $outPath/`printf 'count_%s.txt' $period`
	else if ( $file == kinSIDIS.pdf ) then
	    set fich = $outPath/`printf 'kinSIDIS_%s.pdf' $period`
	else
	    set fich = $outPath/`printf 'kinSIDIS_%s.root' $period`
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
printf '\ * analySIDIS: Script execution status = %d\n' $error
printf '+++++++++++++++++++++++++++++++\n'
exit $copyStatus
