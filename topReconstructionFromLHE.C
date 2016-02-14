#include "topReconstructionFromLHE.h"

using namespace std;

int main()
{
    topReconstructionFromLHE t;
    t.Loop("output_files", 0, 1);

    t.Plot("plots");

    return 0;
}
