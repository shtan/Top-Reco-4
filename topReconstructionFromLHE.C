#include "topReconstructionFromLHE.h"

using namespace std;

int main()
{
    topReconstructionFromLHE t;
    t.debug_verbosity = 2;
    t.Loop("output_files", 0, 1, 10);

//     t.Plot("plots");
    //t.Print();

    return 0;
}
