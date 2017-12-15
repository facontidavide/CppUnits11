#pragma once
#include <ratio>
#include <cmath>

namespace Units{


namespace internal {

struct BaseQuantity{};

// The "RQuantity" class is the prototype template container class, that just holds a double value. The
// class SHOULD NOT BE INSTANTIATED directly by itself, rather use the quantity types defined below.
template<typename MassDim, typename LengthDim, typename TimeDim, typename QAngleDim>
class RQuantity
{
private:
    double value;
    constexpr RQuantity() : value(0.0) {}
    constexpr RQuantity(double val) : value(val) {}

public:

    constexpr RQuantity(double val, BaseQuantity) : value(val) {}

    // The intrinsic operations for a quantity with a unit is addition and subtraction
    constexpr RQuantity const& operator+=(const RQuantity& rhs)
    {
        value += rhs.value;
        return *this;
    }
    constexpr RQuantity const& operator-=(const RQuantity& rhs)
    {
        value -= rhs.value;
        return *this;
    }

    // Returns the value of the quantity in multiples of the specified unit
    constexpr double Convert(const RQuantity& rhs) const
    {
        return value / rhs.value;
    }

    // returns the raw value of the quantity (should not be used)
    constexpr double getValue() const
    {
        return value;
    }
};

template <typename M, typename L, typename T, typename A>
constexpr internal::RQuantity<M, L, T, A>
    operator+(const RQuantity<M, L, T, A>& lhs, const RQuantity<M, L, T, A>& rhs)
{
    return RQuantity<M, L, T, A>(lhs.getValue() + rhs.getValue(), internal::BaseQuantity());
}
template <typename M, typename L, typename T, typename A>
constexpr internal::RQuantity<M, L, T, A>
    operator-(const RQuantity<M, L, T, A>& lhs, const RQuantity<M, L, T, A>& rhs)
{
    return RQuantity<M, L, T, A>(lhs.getValue() - rhs.getValue());
}
template <typename M1, typename L1, typename T1, typename A1,
          typename M2, typename L2, typename T2, typename A2>
constexpr RQuantity<std::ratio_add<M1, M2>, std::ratio_add<L1, L2>,
                    std::ratio_add<T1, T2>, std::ratio_add<A1, A2>>
    operator*(const RQuantity<M1, L1, T1, A1>& lhs, const RQuantity<M2, L2, T2, A2>& rhs)
{
    return RQuantity<std::ratio_add<M1, M2>, std::ratio_add<L1, L2>,
                     std::ratio_add<T1, T2>, std::ratio_add<A1, A2>>
                    (lhs.getValue()*rhs.getValue(), internal::BaseQuantity() );
}
template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A>
    operator*(const double& lhs, const RQuantity<M, L, T, A>& rhs)
{
    return RQuantity<M, L, T, A>(lhs*rhs.getValue(), internal::BaseQuantity() );
}
template <typename M1, typename L1, typename T1, typename A1,
          typename M2, typename L2, typename T2, typename A2>
constexpr RQuantity<std::ratio_subtract<M1, M2>, std::ratio_subtract<L1, L2>,
                    std::ratio_subtract<T1, T2>, std::ratio_subtract<A1, A2>>
    operator/(const RQuantity<M1, L1, T1, A1>& lhs, const RQuantity<M2, L2, T2, A2>& rhs)
{
    return RQuantity<std::ratio_subtract<M1, M2>, std::ratio_subtract<L1, L2>,
                     std::ratio_subtract<T1, T2>, std::ratio_subtract<A1, A2>>
                    (lhs.getValue() / rhs.getValue(), internal::BaseQuantity() );
}
template <typename M, typename L, typename T, typename A>
constexpr RQuantity<std::ratio_subtract<std::ratio<0>, M>, std::ratio_subtract<std::ratio<0>, L>,
                    std::ratio_subtract<std::ratio<0>, T>, std::ratio_subtract<std::ratio<0>, A>>
    operator/(double x, const RQuantity<M, L, T, A>& rhs)
{
    return RQuantity<std::ratio_subtract<std::ratio<0>, M>, std::ratio_subtract<std::ratio<0>, L>,
                     std::ratio_subtract<std::ratio<0>, T>, std::ratio_subtract<std::ratio<0>, A>>
                    (x / rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr RQuantity<M, L, T, A>
    operator/(const RQuantity<M, L, T, A>& rhs, double x)
{
    return RQuantity<M, L, T, A>(rhs.getValue() / x, internal::BaseQuantity() );
}

// Comparison operators for quantities:
// ------------------------------------
template <typename M, typename L, typename T, typename A>
constexpr bool operator==(const RQuantity<M, L, T, A>& lhs, const RQuantity<M, L, T, A>& rhs)
{
    return (lhs.getValue() == rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator!=(const RQuantity<M, L, T, A>& lhs, const RQuantity<M, L, T, A>& rhs)
{
    return (lhs.getValue() != rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator<=(const RQuantity<M, L, T, A>& lhs, const RQuantity<M, L, T, A>& rhs)
{
    return (lhs.getValue() <= rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator>=(const RQuantity<M, L, T, A>& lhs, const RQuantity<M, L, T, A>& rhs)
{
    return (lhs.getValue() >= rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator< (const RQuantity<M, L, T, A>& lhs, const RQuantity<M, L, T, A>& rhs)
{
    return (lhs.getValue()<rhs.getValue());
}
template <typename M, typename L, typename T, typename A>
constexpr bool operator> (const RQuantity<M, L, T, A>& lhs, const RQuantity<M, L, T, A>& rhs)
{
    return (lhs.getValue()>rhs.getValue());
}

} // end namespace internal

// Predefined (physical unit) quantity types:
// ------------------------------------------
#define QUANTITY_TYPE(_Mdim, _Ldim, _Tdim, _Adim, name) \
    typedef internal::RQuantity<std::ratio<_Mdim>, std::ratio<_Ldim>, std::ratio<_Tdim>, std::ratio<_Adim>> name;

// Replacement of "double" type
QUANTITY_TYPE(0, 0, 0, 0, Number);

// Physical quantity types
QUANTITY_TYPE(1, 0, 0, 0, QMass);
QUANTITY_TYPE(0, 1, 0, 0, QLength);
QUANTITY_TYPE(0, 2, 0, 0, QArea);
QUANTITY_TYPE(0, 3, 0, 0, QVolume);
QUANTITY_TYPE(0, 0, 1, 0, QTime);
QUANTITY_TYPE(0, 1, -1, 0, QSpeed);
QUANTITY_TYPE(0, 1, -2, 0, QAcceleration);
QUANTITY_TYPE(0, 1, -3, 0, QJerk);
QUANTITY_TYPE(0, 0, -1, 0, QFrequency);
QUANTITY_TYPE(1, 1, -2, 0, QForce);
QUANTITY_TYPE(1, -1, -2, 0, QPressure);

// QAngle type:
QUANTITY_TYPE(0, 0, 0, 1, QAngle);



// Predefined units:
// -----------------

// Predefined mass units:
constexpr QMass Kg(1.0, internal::BaseQuantity());                            // SI base unit
constexpr QMass Gramme = 0.001 * Kg;
constexpr QMass Tonne = 1000 * Kg;
constexpr QMass Ounce = 0.028349523125 * Kg;
constexpr QMass Pound = 16 * Ounce;
constexpr QMass Stone = 14 * Pound;

// Predefined length-derived units
constexpr QLength Metre(1.0, internal::BaseQuantity());                   // SI base unit
constexpr QLength Decimetre = Metre / 10;
constexpr QLength Centimetre = Metre / 100;
constexpr QLength Millimetre = Metre / 1000;
constexpr QLength Kilometre = 1000 * Metre;
constexpr QLength Inch = 2.54 * Centimetre;
constexpr QLength Foot = 12 * Inch;
constexpr QLength Yard = 3 * Foot;
constexpr QLength Mile = 5280 * Foot;

constexpr QArea Kilometre2 = Kilometre*Kilometre;
constexpr QArea Metre2 = Metre*Metre;
constexpr QArea Decimetre2 = Decimetre*Decimetre;
constexpr QArea Centimetre2 = Centimetre*Centimetre;
constexpr QArea Millimetre2 = Millimetre * Millimetre;
constexpr QArea Inch2 = Inch*Inch;
constexpr QArea Foot2 = Foot*Foot;
constexpr QArea Mile2 = Mile*Mile;

constexpr QVolume Kilometre3 = Kilometre2*Kilometre;
constexpr QVolume Metre3 = Metre2*Metre;
constexpr QVolume Decimetre3 = Decimetre2*Decimetre;
constexpr QVolume Litre = Decimetre3;
constexpr QVolume Centimetre3 = Centimetre2*Centimetre;
constexpr QVolume Millimetre3 = Millimetre2 * Millimetre;
constexpr QVolume Inch3 = Inch2*Inch;
constexpr QVolume Foot3 = Foot2*Foot;
constexpr QVolume Mile3 = Mile2*Mile;

// Predefined time-derived units:
constexpr QTime Second(1.0, internal::BaseQuantity());                        // SI base unit
constexpr QTime Minute = 60 * Second;
constexpr QTime Hour = 60 * Minute;
constexpr QTime Day = 24 * Hour;

constexpr QFrequency Hz(1.0, internal::BaseQuantity());

// Predefined mixed units:
constexpr QAcceleration G = 9.80665 *  Metre / (Second*Second);

constexpr QForce Newton(1.0, internal::BaseQuantity());
constexpr QForce Poundforce = Pound*G;
constexpr QForce Kilopond = Kg*G;

constexpr QPressure Pascal(1.0, internal::BaseQuantity());
constexpr QPressure Bar = 100000 * Pascal;
constexpr QPressure Psi = Pound*G / Inch2;


// Physical unit literals:
// -----------------------

// literals for length units
constexpr QLength operator"" _mm(long double x) { return static_cast<double>(x)*Millimetre; }
constexpr QLength operator"" _cm(long double x) { return static_cast<double>(x)*Centimetre; }
constexpr QLength operator"" _m(long double x) { return static_cast<double>(x)*Metre; }
constexpr QLength operator"" _km(long double x) { return static_cast<double>(x)*Kilometre; }
constexpr QLength operator"" _mi(long double x) { return static_cast<double>(x)*Mile; }
constexpr QLength operator"" _yd(long double x) { return static_cast<double>(x)*Yard; }
constexpr QLength operator"" _ft(long double x) { return static_cast<double>(x)*Foot; }
constexpr QLength operator"" _in(long double x) { return static_cast<double>(x)*Inch; }
constexpr QLength operator"" _mm(unsigned long long int x) { return static_cast<double>(x)*Millimetre; }
constexpr QLength operator"" _cm(unsigned long long int  x) { return static_cast<double>(x)*Centimetre; }
constexpr QLength operator"" _m(unsigned long long int  x) { return static_cast<double>(x)*Metre; }
constexpr QLength operator"" _km(unsigned long long int  x) { return static_cast<double>(x)*Kilometre; }
constexpr QLength operator"" _mi(unsigned long long int  x) { return static_cast<double>(x)*Mile; }
constexpr QLength operator"" _yd(unsigned long long int  x) { return static_cast<double>(x)*Yard; }
constexpr QLength operator"" _ft(unsigned long long int  x) { return static_cast<double>(x)*Foot; }
constexpr QLength operator"" _in(unsigned long long int  x) { return static_cast<double>(x)*Inch; }

// literals for speed units
constexpr QSpeed operator"" _mps(long double x) { return QSpeed(x, internal::BaseQuantity()); }
constexpr QSpeed operator"" _miph(long double x) { return static_cast<double>(x)*Mile / Hour; }
constexpr QSpeed operator"" _kmph(long double x) { return static_cast<double>(x)*Kilometre / Hour; }
constexpr QSpeed operator"" _mps(unsigned long long int x)
                                { return QSpeed(static_cast<long double>(x), internal::BaseQuantity()); }
constexpr QSpeed operator"" _miph(unsigned long long int x)
                                 { return static_cast<double>(x)*Mile / Hour; }
constexpr QSpeed operator"" _kmph(unsigned long long int x)
                                 { return static_cast<double>(x)*Kilometre / Hour; }

// literal for frequency unit
constexpr QFrequency operator"" _Hz(long double x) { return QFrequency(x, internal::BaseQuantity()); }
constexpr QFrequency operator"" _Hz(unsigned long long int x)
                                   { return QFrequency(static_cast<long double>(x), internal::BaseQuantity()); }

// literals for time units
constexpr QTime operator"" _s(long double x) { return QTime(x, internal::BaseQuantity()); }
constexpr QTime operator"" _min(long double x) { return static_cast<double>(x)*Minute; }
constexpr QTime operator"" _h(long double x) { return static_cast<double>(x)*Hour; }
constexpr QTime operator"" _day(long double x) { return static_cast<double>(x)*Day; }
constexpr QTime operator"" _s(unsigned long long int x) { return QTime(static_cast<double>(x), internal::BaseQuantity()); }
constexpr QTime operator"" _min(unsigned long long int x) { return static_cast<double>(x)*Minute; }
constexpr QTime operator"" _h(unsigned long long int x) { return static_cast<double>(x)*Hour; }
constexpr QTime operator"" _day(unsigned long long int x) { return static_cast<double>(x)*Day; }

// literals for mass units
constexpr QMass operator"" _kg(long double x) { return QMass(x, internal::BaseQuantity()); }
constexpr QMass operator"" _g(long double x) { return static_cast<double>(x)*Gramme; }
constexpr QMass operator"" _t(long double x) { return static_cast<double>(x)*Tonne; }
constexpr QMass operator"" _oz(long double x) { return static_cast<double>(x)*Ounce; }
constexpr QMass operator"" _lb(long double x) { return static_cast<double>(x)*Pound; }
constexpr QMass operator"" _st(long double x) { return static_cast<double>(x)*Stone; }
constexpr QMass operator"" _kg(unsigned long long int x) { return QMass(static_cast<double>(x), internal::BaseQuantity()); }
constexpr QMass operator"" _g(unsigned long long int x) { return static_cast<double>(x)*Gramme; }
constexpr QMass operator"" _t(unsigned long long int x) { return static_cast<double>(x)*Tonne; }
constexpr QMass operator"" _oz(unsigned long long int x) { return static_cast<double>(x)*Ounce; }
constexpr QMass operator"" _lb(unsigned long long int x) { return static_cast<double>(x)*Pound; }
constexpr QMass operator"" _st(unsigned long long int x) { return static_cast<double>(x)*Stone; }

// literals for acceleration units
constexpr QAcceleration operator"" _mps2(long double x) { return QAcceleration(x, internal::BaseQuantity()); }
constexpr QAcceleration operator"" _mps2(unsigned long long int x)
                                        { return QAcceleration(static_cast<double>(x), internal::BaseQuantity()); }
constexpr QAcceleration operator"" _G(long double x) { return static_cast<double>(x)*G; }
constexpr QAcceleration operator"" _G(unsigned long long int x) { return static_cast<double>(x)*G; }

// literals for force units
constexpr QForce operator"" _N(long double x) { return QForce(x, internal::BaseQuantity()); }
constexpr QForce operator"" _N(unsigned long long int x) { return QForce(static_cast<double>(x), internal::BaseQuantity()); }
constexpr QForce operator"" _lbf(long double x) { return static_cast<double>(x)*Poundforce; }
constexpr QForce operator"" _lbf(unsigned long long int x) { return static_cast<double>(x)*Poundforce; }
constexpr QForce operator"" _kp(long double x) { return static_cast<double>(x)*Kilopond; }
constexpr QForce operator"" _kp(unsigned long long int x) { return static_cast<double>(x)*Kilopond; }

// literals for pressure units
constexpr QPressure operator"" _Pa(long double x) { return QPressure(x, internal::BaseQuantity()); }
constexpr QPressure operator"" _Pa(unsigned long long int x)
                                  { return QPressure(static_cast<double>(x), internal::BaseQuantity()); }
constexpr QPressure operator"" _bar(long double x) { return static_cast<double>(x)*Bar; }
constexpr QPressure operator"" _bar(unsigned long long int x) { return static_cast<double>(x)*Bar; }
constexpr QPressure operator"" _psi(long double x) { return static_cast<double>(x)*Psi; }
constexpr QPressure operator"" _psi(unsigned long long int x) { return static_cast<double>(x)*Psi; }


// Angular unit literals:
// ----------------------
constexpr QAngle Radian(1.0, internal::BaseQuantity());

constexpr QAngle operator"" _pi(long double x)
    { return static_cast<double>(x) * 3.1415926535897932384626433832795 * Radian; }
constexpr QAngle operator"" _pi(unsigned long long int x)
    { return (static_cast<double>(x) * 3.1415926535897932384626433832795) * Radian; }

// Predefined angle units:
constexpr QAngle Degree = (2.0_pi / 360.0);

// literals for angle units
constexpr QAngle operator"" _rad(long double x) { return QAngle(x, internal::BaseQuantity()); }
constexpr QAngle operator"" _rad(unsigned long long int x) { return QAngle(static_cast<double>(x), internal::BaseQuantity()); }
constexpr QAngle operator"" _deg(long double x) { return static_cast<double>(x)*Degree; }
constexpr QAngle operator"" _deg(unsigned long long int x) { return static_cast<double>(x)*Degree; }

// Conversion macro, which utilizes the string literals
#define ConvertTo(_x, _y) (_x).Convert(1.0_##_y)



// Typesafe mathematical operations:
// ---------------------------------
template <typename M, typename L, typename T, typename A>
constexpr internal::RQuantity<std::ratio_divide<M, std::ratio<2>>, std::ratio_divide<L, std::ratio<2>>,
                    std::ratio_divide<T, std::ratio<2>>, std::ratio_divide<A, std::ratio<2>>>
    Qsqrt(const internal::RQuantity<M, L, T, A>& num)
{
    return internal::RQuantity<std::ratio_divide<M, std::ratio<2>>, std::ratio_divide<L, std::ratio<2>>,
                     std::ratio_divide<T, std::ratio<2>>, std::ratio_divide<A, std::ratio<2>>>
                    (sqrt(num.getValue()));
}

// Typesafe trigonometric operations
inline double sin(const QAngle &num)
{
    return std::sin(num.getValue());
}
inline double cos(const QAngle &num)
{
    return std::cos(num.getValue());
}
inline double tan(const QAngle &num)
{
    return std::tan(num.getValue());
}

} // end namespace Units


