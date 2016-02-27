#include "topReconstructionFromLHE.h"

using namespace std;

int main()
{
    topReconstructionFromLHE t;
    t.debug_verbosity = 1;
    t.Loop("output_files", 0, 1, 5);

//     t.Plot("plots");

    return 0;
}
