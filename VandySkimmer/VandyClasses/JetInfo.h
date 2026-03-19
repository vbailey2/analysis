#include <iostream>
#include <limits>
#include <vector>

class JetInfo {
    public:
        JetInfo() = default;
        ~JetInfo();

        double px() {return Px;};
        void set_px(double tpx) {Px = tpx;};
        
        double py() {return Py;};
        void set_py(double tpy) {Py = tpy;};

        double pz() {return Pz;};
        void set_pz(double tpz) {Pz = tpz;};

        double e() {return E;};
        void set_e(double te) {E = te;};

        double pt() {return Pt;};
        void set_pt(double pt) {Pt = pt;};

        double emCaloFrac() {return EMCaloFrac;};
        void set_emCaloFrac(double cf) {EMCaloFrac = cf;};

        std::vector<int> get_constituents() { return constituents; };
        void set_constituents (std::vector<int> jetCons) { constituents = jetCons; };

        void CopyTo(JetInfo *jet);
    private:
        double Px;
        double Py;
        double Pz;
        double E;
        double Pt;
        double EMCaloFrac;
        std::vector<int> constituents;
};
