#include "topReconstructionFromLHE.h"

using namespace std;
//using namespace commonstruct;

int main()
{
    topReconstructionFromLHE t;
    t.debug_verbosity = 1;
    //t.Loop("output_files", 0, 1, 1, 10);
    
    //t.Loop("testoutput", 0, 10, 0, 10);

//     t.Plot("plots");
    t.Print();
/*
    top_system foxtrot(26);
    //foxtrot.input.b_px = 5;
    foxtrot.vars.b_delta_pt = 8;
    cout << foxtrot.vars.b_pt() << endl;

    foxtrot.input.btest = 6.7;
    cout << foxtrot.vars.b_pt() << endl;

    hadronic_top_system waltz(134, 5.3);
    waltz.vars.b_delta_pt = 2;
    cout << waltz.vars.b_pt() << endl;
    cout << waltz.optional1_ << endl;

    cout << foxtrot.best.lars() << endl;

    cout << adder(foxtrot) << endl;
    foxtrot.best.bestchi = adder(foxtrot);
    cout << foxtrot.best.bestchi << endl;
*/
    return 0;
}
