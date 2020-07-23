#ifndef __COMPMULT_H_INCLUDED__
#define __COMPMULT_H_INCLUDED__

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <TLatex.h>
#include <TLegend.h>

using namespace std;

struct Acceptance { double tab[2][2][2][4]; };

Acceptance fAcceptance1[9][6][12];
Acceptance fAcceptance2[9][6][12];
Acceptance fAcceptance1_yavg[9][12];
Acceptance fAcceptance2_yavg[9][12];
Acceptance fAcceptance1_zavg[9];
Acceptance fAcceptance2_zavg[9];

double z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};

//Graphic Style

int fMarkerColor[6] = {2,95,209,226,4,221};
int fMarkerStyle[6][2] = {{24,20},{26,22},{25,21},{27,33},{28,34},{30,29}};

#endif
