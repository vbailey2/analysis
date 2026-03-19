#include "JetInfo.h"

JetInfo::~JetInfo() = default;

void JetInfo::CopyTo(JetInfo* jet)
{
    jet->set_px(Px);
    jet->set_py(Py);
    jet->set_pz(Pz);
    jet->set_e(E);
    jet->set_pt(Pt);
    jet->set_emCaloFrac(EMCaloFrac);
    std::vector<int> constit;
    for(auto cons : constituents)
    {
        constit.push_back(cons);
    }
    jet->set_constituents(constit);

    return;
}
