# Copyright (C) 2012-2018 The ESPResSo project
# Copyright (C) 2012 Olaf Lenz
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

git ls-files --exclude-standard |
egrep -v '\.(blk|gz|data|dat|tab|chk|jpg|png|pdf|fig|gif|xcf|bib|vtf|vtk|svg|ico|eps)$' |
egrep -v '^testsuite/configs/|^old/|^cmake/|^libs/' |
egrep -v '(ChangeLog|AUTHORS|COPYING|NEWS|README|INSTALL)' |
egrep -v '(\.gitignore|pkgIndex\.tcl)' |
egrep -v '(config/config\.guess|config/config\.sub|config/install-sh)' |
egrep -v '(Doxyfile|latexmk\.1|latexmkrc|assemble_quickref\.awk|doc/misc/homepage/palette\.html)' |
egrep -v '(src/features\.def)' |
egrep -v '(featurelist)' |
egrep -v '(\.cproject|\.project|\.settings)'
