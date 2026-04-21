// Copyright (C) 2026 Wenbo Lyu
//
// This file is part of pc.
//
// pc is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// pc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with infoxtr. If not, see <https://www.gnu.org/licenses/>.

#ifndef PC_PC_H
#define PC_PC_H

// ============================================================
// Dependency Guard: Encourage best practices (non-blocking)
// ============================================================

#if defined(Rcpp_hpp) && !defined(COMPILING_PC)
    #warning "It is recommended to include <pc.h> alone, as it already includes <Rcpp.h> and <RcppThread.h>."
#endif

// ============================================================
// Core Dependencies (Auto-included for users)
// ============================================================

#include <Rcpp.h>
#include <RcppThread.h>

// ============================================================
// Module Headers (Organized by functionality)
// ============================================================

#include "pc/embed.hpp"
#include "pc/combn.hpp"
#include "pc/numericutils.hpp"
#include "pc/distance.hpp"
#include "pc/symdync.hpp"

// ============================================================
// Convenience Converters (Inline helpers for R/C++ interop)
// ============================================================

#endif // PC_PC_H
