/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANB WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Barycentric<Cmpt>::Barycentric
(
    const Cmpt& va,
    const Cmpt& vb,
    const Cmpt& vc,
    const Cmpt& vd
)
{
    this->v_[A] = va;
    this->v_[B] = vb;
    this->v_[C] = vc;
    this->v_[D] = vd;
}


//template<class Cmpt>
//inline Foam::Barycentric<Cmpt>::Barycentric
//(
//    const Cmpt& va,
//    const Cmpt& vb,
//    const Cmpt& vc
//)
//:
//    Barycentric<Cmpt>(va, vb, vc, 1 - va - vb - vc)
//{
//}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::Barycentric<Cmpt>::a() const
{
    return this->v_[A];
}


template<class Cmpt>
inline const Cmpt& Foam::Barycentric<Cmpt>::b() const
{
    return this->v_[B];
}


template<class Cmpt>
inline const Cmpt& Foam::Barycentric<Cmpt>::c() const
{
    return this->v_[C];
}


template<class Cmpt>
inline const Cmpt& Foam::Barycentric<Cmpt>::d() const
{
    return this->v_[D];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric<Cmpt>::a()
{
    return this->v_[A];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric<Cmpt>::b()
{
    return this->v_[B];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric<Cmpt>::c()
{
    return this->v_[C];
}


template<class Cmpt>
inline Cmpt& Foam::Barycentric<Cmpt>::d()
{
    return this->v_[D];
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
inline Barycentric<Cmpt> toBarycentric
(
    const Tensor<Cmpt>& invA,
    const Vector<Cmpt>& o,
    const Vector<Cmpt>& v
)
{
    const Vector<Cmpt> w = invA&(v - o);
    return Barycentric<Cmpt>(w.x(), w.y(), w.z(), 1 - w.x() - w.y() - w.z());
}


template<class Cmpt>
inline Barycentric<Cmpt> toBarycentric
(
    const Tensor<Cmpt>& invA,
    const Vector<Cmpt>& v
)
{
    const Vector<Cmpt> w = invA&v;
    return Barycentric<Cmpt>(w.x(), w.y(), w.z(), - w.x() - w.y() - w.z());
}


template<class Cmpt>
inline Vector<Cmpt> fromBarycentric
(
    const Tensor<Cmpt>& A,
    const Vector<Cmpt>& o,
    const Barycentric<Cmpt>& b
)
{
    const Vector<Cmpt> v(b.a(), b.b(), b.c());
    return (A&v) + o;
}


template<class Cmpt>
inline Vector<Cmpt> fromBarycentric
(
    const Tensor<Cmpt>& A,
    const Barycentric<Cmpt>& b
)
{
    const Vector<Cmpt> v(b.a(), b.b(), b.c());
    return A&v;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
