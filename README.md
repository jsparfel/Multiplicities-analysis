This is a copy of "github.com/nipierre/Multiplicities.git".

# Multiplicities Analysis

## Summary
 1. [Building](#Building)
 2. [Flow](#Flow)
 3. [analySIDIS_split](#analySIDIS-split)
 4. [analySIDIS_collect](#analySIDIS-collect)
 5. [acceptance_split](#acceptance-split)
 6. [acceptance_fuse](#acceptance-fuse)
 7. [acceptance_collect](#acceptance-collect)
 8. [submitJobs.py](#submitJobs.py)
 9. [compMult](#compMult)
 10. [compAcc](#compAcc)

## Building

 MAKEFILE PHONY TARGETS:
  - `make`: builds all targets
  - `make analySIDIS`: builds the SIDIS analysis framework.
  - `make acceptance`: builds the acceptance analysis framework
  - `make comparison`: builds several comparison applets
  - `make clean`

## Flow

To calculate raw multiplicities: **Multiplicities Calculation** _[analySIDIS_split --> analySIDIS_collect]_

To calculate acceptance and corrected multiplicities: **Acceptance Calculation** _[acceptance_split --> acceptance_fuse --> acceptance_collect]_ **-->** **Multiplicities Calculation** _[analySIDIS_split --> analySIDIS_collect]_


## [C++] analySIDIS_split<a name="analySIDIS-split" />

**Description:**
Takes the TTree and does the cut of the analysis. Outputs DIS Event and Hadron counts.

**Requires:**
 - **ROOT TTree from Phast User Event 120** 
 - **RICH matrices in `data/rich_mat_2016_julien.txt` and `data/rich_mat_error_julien.txt`**

**User Dependence in file (see beginning of `analySIDIS_split.cc`):**
 - **data_path to Real Data: e.g. `#define data_path "/sps/compass/julien"`**
 - **RICH matrices to be used: e.g. `#define mat_RICH_2016_name "data/rich_mat_2016_julien.txt"` and `#define err_RICH_name "data/rich_mat_error_julien.txt"`**

**In File Flags (see beginning of `analySIDIS_split.cc`):**
 - **Momentum Boundaries**
 - **X Boundaries**
 - **Y Boundaries**
 - **Muon charge separation [YES/NO | 1/0]**
 - **RICH matrix correction or not: e.g. `#define RICH 1`**

**Directories/files to be created before execution:**
 - **`[data_path]/P07/`, `[data_path]/P08/` etc. where `P07`, `P08`, etc. are directories containing .root files of the corresponding period, previously treated with Phast User Event 120 and `[data_path]` is defined in `analySIDIS_split.cc` (see above)**
 - **`[data_path]/P07/filelist.txt` where `filelist.txt` is a text file containing the list of the names of .root files to be treated for the corresponding period**
 - **`./PeriodFile.txt` where `.` denotes the directory from which the execution of `analySIDIS_split` is launched and `PeriodFile.txt` is a text file with structure defined below**
 - **`./data/` contains files such as `rich_mat_2016_julien.txt` etc.**
 - **`./rawmult/2016/` is the directory where output files will be created**
 

**Call:**
```Bash
[exe_path]/analySIDIS_split [PeriodFile] [OPTIONAL FLAG]
```
where `[exe_path]` is a directory containing the executable `analySIDIS_split` and `[PeriodFile]` is a file containing the periods to treat. It has the following structure:

```
P01 0
P02 0
P03 0
P04 0
P05 0
P06 0
P07 1
P08 1
P09 1
P10 1
P11 1
```

`[OPTIONAL FLAG]`
 - `-k` : draw kinematic plots

See `submitJobs.py` to submit `analySIDIS_split` jobs to the batch (necessary when treating massive data).

**Outputs:**
 - **Count files in `./rawmult/2016/`**
 - **If option `-k` is passed, kinematic plots `kinSIDIS.pdf` and `kinSIDIS.root` in `./`**

## [C++] analySIDIS_collect<a name="analySIDIS-collect" />

**Description:**
Takes the output of `analySIDIS_split` (and of `acceptance_collect` if acceptance correction is activated), computes the Multiplicities and stores/plots them.

**Requires:**
 - **Output of `analySIDIS_split` in `./rawmult/2016/` where `.` denotes the directory from which `analySIDIS_collect` is launched**
 - **If acceptance correction activated, output from `acceptance_collect` in `./Multiplicities/acceptance/2016/`**
 - **Inclusive and Semi-inclusive Radiative Corrections: file `./data/proton_semi_inclusive_RC.txt`**
 - **Diffractive Vector Meson Correction: file `./data/DVM_2016.dat`**

**In File Flags**
 - **Diffractive Vector Meson Correction [YES/NO | 1/0] e.g. `#define DVMC 0`**
 - **Radiative Corrections [YES/NO | 1/0] e.g. `#define SIRC 0`**
 - **No Acceptance [YES/NO | 1/0] e.g. `#define NO_ACC 1` to de-activate acceptance correction (raw multiplicities)**
 - **Y Integration or mean [Mean/Weighted Mean/Integration | 1/2/3] e.g. `#define YMULT 2` (my advice is to use 2)**
 - **Staggered Multiplicities in plot [YES/NO | 1/0] e.g. `#define STAGGERED 1` (if set to 0 multiplicities in different bins of y will not be staggered and hence may not be distinguishable)**
 - **Data path in which output files will be created: e.g. `#define data_path "./Multiplicities"`**

**Directories/files to be created before execution:**
 - **`./PeriodFile.txt` is a text file with structure defined above**
 - **`./data/` is a directory containing files such as `proton_semi_inclusive_RC.txt` etc.**
 - **`./Multiplicities/` is the directory where output files will be created**

**Call:**
 ```Bash
 [exe_path]/analySIDIS_collect [PeriodFile]
 ```
where `[exe_path]` is a directory containing the executable `analySIDIS_collect` and `[PeriodFile]` is a file containing the periods to treat (see above in analySIDIS_split for the structure of this file).

**Outputs:**
 - **Multiplicity text files in `./Multiplicities`:**
  - **`multiplicities_{}.txt` for Hadron, Pion, Kaon and Proton**
  - **`multiplicities_{}_yavg.txt` for Hadron, Pion, Kaon and Proton**
  - **`multiplicities_raw.txt`**
  - **`multiplicities_h{}.txt` for + and -**
  - **`multiplicities_hadron_{}.txt` for pt and theta**
  - **`reldiff.txt`**
 - **Multiplicity plots in `data_path`:**
  - **`{}_multiplicity_file.pdf` for Hadron, Pion, Kaon and Proton: plots multiplicities in bins of x,y,z**
  - **`{}_multiplicity_zvtx_file.pdf` for Hadron, Pion, Kaon and Proton: plots multiplicities for the different bins of Zvtx, allowing to check stability over target**
  - **`{}_multiplicity_period_file.pdf` for Hadron, Pion, Kaon and Proton: plots multiplicities for the different periods, allowing to check stability over time**
  - **`{}_multiplicity_yavg_file.pdf` for Hadron, Pion, Kaon and Proton: plots y-averaged multiplicities**
  - **`{}_multiplicity_sum_file.pdf` for Hadron, Pion, Kaon and Proton: plots y-averaged and z-integrated multiplicities sum**
  - **`{}_multiplicity_ratio_file.pdf` for Hadron, Pion, Kaon and Proton: plots y-averaged and z-integrated multiplicities ratio**

## [C++] acceptance_split<a name="acceptance-split" />


**Description:**
Takes the Monte Carlo .root files and does the cut of the analysis. Outputs DIS event and Hadron counts. Important: acceptance_split MUST be launched separately for mu+ and mu- .root MC files and the outputs of these 2 executions should be put in 2 different directories (for the output use a symbolic link `acceptance` pointing alternatively towards a directory for mu+ and for mu-).

**User Dependence in file (see beginning of acceptance_split.cc)**
 - **data_path to MC data: e.g. `#define data_path "/sps/compass/julien/MC"`**

**In File Flags:**
 - **Data year: e.g. `#define Y2016 1` to treat 2016 data**
 - **Momentum Boundaries**
 - **X Boundaries**
 - **Y Boundaries**
 - **W Boundaries**
 - **XX0 Limit**

**Directories/files to be created before execution:**
 - **`[data_path]/P07/`, `[data_path]/P08/` etc. where `P07`, `P08`, etc. are directories containing .root files for MC data of the corresponding period, previously treated with Phast User Event 120 and `[data_path]` is defined in `acceptance_split.cc` (see above)**
 - **`[data_path]/P07/filelist.txt` where `filelist.txt` is a text file containing the list of the names of .root files to be treated for the corresponding period**
 - **`./PeriodFile.txt` where `.` denotes the directory from which the execution of `accsplit` is launched and `PeriodFile.txt` is a text file with structure defined above**
 - **`./acceptance/2016/DIS/`, `./acceptance/2016/electron/` and `./acceptance/2016/hadron/` are directories into which output files will be created**

**Call:**
```Bash
[exe_path]/accsplit [PeriodFile] [OPTIONAL FLAG]
```
where `[exe_path]` is a directory containing the executable `accsplit` and `[PeriodFile]` is a file containing the periods to treat.

`[OPTIONAL FLAG]`
 - `-k` : draw kinematic plots

See `submitJobs.py` to submit `acceptance_split` jobs to the batch (necessary when treating massive data).

**Outputs:**
 - **Count files in `./acceptance/2016/`**
 - **If option `-k` is passed, kinematic plots `kinMC.pdf` and `kinMC.root` in `./`**

## [C++] acceptance_fuse<a name="acceptance-fuse" />


**Description:**
Fuses mu+/- counts from the 2 different executions of acceptance_split (see above).

**Requires:**
 - **Outputs from 2 different executions of `acceptance_split` (one for mu+ .root files and one for mu- .root files)**

**Directories/files to be created before execution:**
 - **`./acceptance/2016/DIS/`, `./acceptance/2016/electron/` and `./acceptance/2016/hadron/` are directories into which output files will be created(again, use a symbolic link `acceptance`)**

**Call:**
```Bash
[exe_path]/accfuse [PeriodName] [MU+ FILELIST] [MU- FILELIST]
```
where `[exe_path]` is a directory containing the executable `accfuse`, `[PeriodName]` is the name of the period, eg. `P09`, and `[MU+ FILELIST]`(`[MU- FILELIST]`) is a file containing the name of the output directory of the mu+(mu-) execution of `acceptance_split`. For example `[MU+ FILELIST]` may contain `accS51Dj+` with the outputs of `acceptance_split` for mu+ being in `./accS51Dj+/2016/`.

**Outputs:**
 - **Fused count files in `./acceptance/2016/`**

## [C++] acceptance_collect<a name="acceptance-collect" />


**Description:**
Takes the output of `acceptance_fuse`, computes the acceptances and stores/plots them.

**Requires:**
 - **Output of `acceptance_fuse` in `./acceptance/2016/` where . is the directory from which `acccollect` is launched**

**User Dependence:**
- **Path to outputs in `acceptance_collect.cc`, e.g. `#define dirroot "/sps/compass/julien/nicolas/Multiplicities/acceptance"`**

**In File Flags**
 - **Data year: e.g. `#define Y2016 1` to treat 2016 data**
 - **Y STAGGERING (SPREAD) [YES/NO | 1/0]**

**Directories/files to be created before execution:**
 - **`./Multiplicities/acceptance/2016/` is the directory into which output files will be created**
 - **`./PeriodFile.txt` is a text file with structure defined above**

**Call:**
```Bash
[exe_path]/acccollect [PeriodFile]
```
where `[exe_path]` is a directory containing the executable `acccollect` and `[PeriodFile]` is a file containing the periods to treat.

**Outputs:**
 - **Acceptance text files in `$dirroot/2016/`:**
  - **`acceptance_{}.txt` per Period**
  - **`acceptance_yavg_{}.txt` per Period**
  - **`acceptance_vtx_{}.txt` per Period**
  - **`acceptance_theta_{}.txt` per Period**
  - **`reldiff_vtx_{}.txt` per Period**
 - **Acceptance plots in `$dirroot/2016/`:**
  - **`{}_acceptance_{}.pdf` for Hadron, Pion, Kaon, Proton and per Period**
  - **`{}_acceptance_corr_{}.pdf` for Hadron, Pion, Kaon, Proton and per Period**

## [Python] submitJobs.py<a name="submitJobs.py" />


**Description:**
Submit `analySIDIS_split` or `acceptance_split` jobs to the batch (one job per period).

**Requires:**
 - **Same as `analySIDIS_split`/`acceptance_split`**
 - **Scripts `analySIDIS.csh` and `acceptance.csh`**

**User Dependence:**
 - **Same as `analySIDIS_split`/`acceptance_split`**

**In File Flags**
 - **Same as `analySIDIS_split`/`acceptance_split`**

**Directories/files to be created before execution:**
 - **Same as `analySIDIS_split`/`acceptance_split`**
 - **Scripts `analySIDIS.csh` and `acceptance.csh` should be put in the same directory as executables `analySIDIS_split` and `accsplit`**

**Call:**
```Bash
python [path]/submitJobs.py -e [exe_path] -s [script_name] -o [output_path] -f [PeriodFile]
```
where `[path]` is a directory containing `submitJobs.py`, `[exe_path]` is a directory containing the executable `analySIDIS_split` or `accsplit`, `[script_name]` is either `acceptance.csh` or `analySIDIS.csh`, `[output_path]` is the directory into which output files will be written and `[PeriodFile]` is a file containing the periods to treat. Option `-f` forces the writing of outputs even if files already exist.

**Outputs:**
 - **Same as `analySIDIS_split`/`acceptance_split`**
 - **`OUTPUT_DIR/P0x.log`: one .log file per job**

## [C++] compMult<a name="compMult" />

**Description:**
Takes two different multiplicities results and plots the ratio in order to compare them. Allows to see if two different calculations of multiplicities are compatible.

**Requires:**
 - **Outputs of 2 different executions of `analySIDIS_collect` in different directories**

**Directories/files to be created before execution:**
 - **`./compMult/` is the directory into which output files will be created**

**Call:**
```Bash
[exe_path]/compMult [MULT_DIRECTORY_1] [MULT_DIRECTORY_2]
```
where `[exe_path]` is a directory containing the executable `compMult` and `[MULT_DIRECTORY_1]`(`[MULT_DIRECTORY_2]`) is the directory where the first(the second) multiplicities files are stored (e.g. `Multiplicities_S51Dj`).

**Outputs:**
 - **Comparison pdf files in `./compMult/`:**
  - **`mult_ratio.pdf` ratio of multiplicities in bins of x,y,z**
  - **`mult_ratio_yavg.pdf` ratio of y-averaged multiplicities in bins of x,z**

## [C++] compAcc<a name="compAcc" />

**Description:**
Takes two different acceptances results and plots the ratio in order to compare them. Allows to see if two different calculations of acceptances are compatible.

**Requires:**
 - **Outputs of 2 different executions of `acceptance_collect` in different directories**

**Directories/files to be created before execution:**
 - **`./compAcc/` is the directory into which output files will be created**

**Call:**
```Bash
[exe_path]/compAcc [ACC_DIRECTORY_1] [ACC_DIRECTORY_2] [PERIOD]
```
where `[exe_path]` is a directory containing the executable `compAcc`, `[ACC_DIRECTORY_1]`(`[ACC_DIRECTORY_2]`) is the directory where the first(the second) acceptances files are stored (e.g. `Multiplicities_S51Dj/acceptance/2016`) and `[PERIOD]` is the name of the period of MC data e.g. `P09`.

**Outputs:**
 - **Comparison pdf files in `./compAcc/`:**
  - **`acc_ratio.pdf` ratio of acceptances in bins of x,y,z**
  - **`acc_ratio_yavg.pdf` ratio of y-averaged acceptances in bins of x,z**
