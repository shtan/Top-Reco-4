#include "topReconstructionFromLHE.h"

using namespace std;

int main()
{
    topReconstructionFromLHE t;
    t.debug = true;
    t.Loop("output_files", 0, 1, 1);

    t.Plot("plots");

    return 0;
}
