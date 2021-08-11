template <class T>
T lj_potential(T const& rij2)
{
    T rm6 = 1./(rij2*rij2*rij2);
    T vlj = (rm6 - 1.0)*rm6;
    return vlj;
}


template <class T>
T lj_force_factor(T const & rij2)
{
    T rm2 = 1.0/rij2;
    T rm6 = (rm2*rm2*rm2);
    T f = (1.0 - 2.0*rm6 )*rm6*rm2*6.0;
    return f;
}
