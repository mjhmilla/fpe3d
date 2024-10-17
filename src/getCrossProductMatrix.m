%%
% SPDX-FileCopyrightText: 2024 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
function vecX = getCrossProductMatrix(vec)

vecX = [0       , -vec(3,1),  vec(2,1);...
        vec(3,1),  0       , -vec(1,1);...
       -vec(2,1),  vec(1,1),  0];