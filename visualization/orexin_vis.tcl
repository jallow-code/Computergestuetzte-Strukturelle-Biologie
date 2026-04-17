# Orexin receptor visualization for VMD
# Structures:
#   6TOS chain A = OX1R + ligand NRE
#   6TPN chain A = OX2R + ligand NU8
#
# Usage in VMD:
#   source orexin_vis.tcl
#   orexin_load
#   scene_ox1
#   scene_ox2
#   scene_overlay
#
# Then render manually from the VMD GUI or with:
#   render TachyonInternal ox1.tga

# Pocket residues often discussed for subtype selectivity
set OX1_POS_2_61 103
set OX1_POS_3_33 127
set OX2_POS_2_61 111
set OX2_POS_3_33 135

set OX1_PDB "6TOS.pdb"
set OX2_PDB "6TPN.pdb"
set OX1_CHAIN "A"
set OX2_CHAIN "A"
set OX1_LIG "NRE"
set OX2_LIG "NU8"

# Explicit OX1R->OX2R CA mapping for the resolved receptor core.
# This avoids fit errors from different loop/fusion coverage.
set core_pairs {
  {27 32} {28 33} {29 35} {30 37} {31 39} {32 40} {33 41} {34 42} {35 43}
  {36 44} {37 45} {38 46} {39 47} {40 48} {41 49} {42 50} {43 51} {44 52}
  {45 53} {46 54} {47 55} {48 56} {49 57} {50 58} {51 59} {52 60} {53 61}
  {54 62} {55 63} {56 64} {57 65} {58 66} {59 67} {60 68} {61 69} {62 70}
  {63 71} {64 72} {65 73} {66 74} {67 75} {68 76} {69 77} {70 78} {71 79}
  {72 80} {73 81} {74 82} {75 83} {76 84} {77 85} {78 86} {79 87} {80 88}
  {81 89} {82 90} {83 91} {84 92} {85 93} {86 94} {87 95} {88 96} {89 97}
  {90 98} {91 99} {92 100} {93 101} {94 102} {95 103} {96 104} {98 105}
  {99 107} {100 108} {101 109} {102 110} {103 111} {104 112} {105 113}
  {106 114} {107 115} {108 116} {109 117} {110 118} {111 119} {112 120}
  {113 121} {114 122} {115 123} {116 124} {117 125} {118 126} {119 127}
  {120 128} {121 129} {122 130} {123 131} {124 132} {125 133} {126 134}
  {127 135} {128 136} {129 137} {130 138} {131 139} {132 140} {133 142}
  {135 143} {136 144} {137 145} {138 146} {139 147} {140 148} {141 149}
  {142 150} {143 151} {144 152} {145 153} {146 154} {147 155} {148 156}
  {149 157} {150 158} {151 159} {152 160} {153 161} {154 162} {155 163}
  {156 164} {157 165} {158 166} {159 167} {160 168} {161 169} {162 170}
  {163 171} {164 172} {165 173} {166 174} {167 175} {168 176} {169 177}
  {170 178} {171 179} {172 180} {173 181} {174 182} {175 183} {176 184}
  {177 185} {178 186} {179 187} {180 188} {181 189} {182 190} {183 191}
  {184 192} {185 193} {186 194} {187 195} {188 196} {197 197} {198 206}
  {199 207} {200 208} {201 209} {202 210} {203 211} {204 212} {205 213}
  {206 214} {207 215} {208 216} {209 217} {210 218} {211 219} {212 220}
  {213 221} {214 222} {215 223} {216 224} {217 225} {218 226} {219 227}
  {220 228} {221 229} {222 230} {223 231} {224 232} {225 233} {226 234}
  {227 235} {228 236} {229 237} {230 238} {231 239} {232 240} {233 241}
  {234 242} {235 243} {236 244} {237 245} {238 246} {239 247} {240 248}
  {241 249} {242 250} {243 251} {244 252} {286 253} {287 254} {288 294}
  {289 295} {290 296} {291 297} {292 298} {293 299} {294 300} {295 301}
  {296 302} {297 303} {298 304} {299 305} {300 306} {301 307} {302 308}
  {303 309} {304 310} {305 311} {306 312} {307 313} {308 314} {309 315}
  {310 316} {311 317} {312 318} {313 319} {314 320} {315 321} {316 322}
  {317 323} {318 324} {319 325} {320 326} {321 327} {322 328} {323 329}
  {324 330} {325 331} {326 332} {327 333} {328 334} {329 335} {330 336}
  {331 337} {332 338} {333 339} {334 340} {335 341} {336 342} {337 343}
  {338 344} {339 345} {340 346} {341 347} {342 348} {343 349} {344 350}
  {345 351} {346 352} {347 353} {348 354} {349 355} {350 356} {351 357}
  {352 358} {353 359} {354 360} {355 361} {356 362} {357 363} {358 364}
  {359 365} {360 366} {361 367} {362 368} {363 369} {364 370} {365 371}
  {366 372} {367 373} {368 374} {369 375} {370 376} {371 377} {372 378}
  {373 379} {374 380} {375 381} {376 382} {377 384}
}

proc _paired_ca_selections {mol1 chain1 mol2 chain2 pairs} {
  set idx1 {}
  set idx2 {}
  foreach pair $pairs {
    set r1 [lindex $pair 0]
    set r2 [lindex $pair 1]

    set s1 [atomselect $mol1 "protein and chain $chain1 and resid $r1 and name CA"]
    set s2 [atomselect $mol2 "protein and chain $chain2 and resid $r2 and name CA"]

    if {[$s1 num] > 0 && [$s2 num] > 0} {
      lappend idx1 [lindex [$s1 get index] 0]
      lappend idx2 [lindex [$s2 get index] 0]
    }

    $s1 delete
    $s2 delete
  }

  if {[llength $idx1] == 0 || [llength $idx2] == 0} {
    error "No matched CA atoms found for alignment. Check chain IDs and residue numbering."
  }

  set sel1 [atomselect $mol1 "index [join $idx1 { }]"]
  set sel2 [atomselect $mol2 "index [join $idx2 { }]"]
  return [list $sel1 $sel2]
}

proc _clear_reps {molid} {
  set nreps [molinfo $molid get numreps]
  for {set i 0} {$i < $nreps} {incr i} {
    mol delrep 0 $molid
  }
}

proc _common_display {} {
  display projection Orthographic
  display depthcue off
  axes location Off
  color Display Background white
  display shadows on
  display ambientocclusion on
  display aoambient 0.85
  display aodirect 0.30
  display nearclip set 0.01
  display farclip set 10.0
  display cuedensity 0.0
  material change ambient AOChalky 0.25
  material change diffuse AOChalky 0.75
  material change specular AOChalky 0.10
  material change shininess AOChalky 0.10
  material change opacity AOChalky 1.00
  material change ambient Transparent 0.25
  material change diffuse Transparent 0.65
  material change specular Transparent 0.05
}

proc _reset_labels {} {
  foreach item [label list Atoms] {
    catch {label delete Atoms $item}
  }
}

proc _add_site_labels {molid chain resid_list} {
  foreach resid $resid_list {
    set sel [atomselect $molid "protein and chain $chain and resid $resid and name CA"]
    if {[$sel num] > 0} {
      label add Atoms $molid/[lindex [$sel get index] 0]
    }
    $sel delete
  }
}

proc _water_shell_selection {chain ligresn cutoff} {
  return "water and same residue as within $cutoff of (chain $chain and resname $ligresn)"
}

proc _base_receptor_selection {chain} {
  return "protein and chain $chain and not resid 1000 to 1199"
}

proc orexin_load {} {
  global OX1_PDB OX2_PDB OX1_CHAIN OX2_CHAIN core_pairs

  mol delete all
  set ox1 [mol new $OX1_PDB type pdb waitfor all]
  set ox2 [mol new $OX2_PDB type pdb waitfor all]

  lassign [_paired_ca_selections $ox1 $OX1_CHAIN $ox2 $OX2_CHAIN $core_pairs] fit1 fit2
  set ox2_all [atomselect $ox2 "all"]
  set trans [measure fit $fit2 $fit1]
  $ox2_all move $trans

  puts "Aligned OX2R onto OX1R using [llength $core_pairs] matched CA positions."
  puts "Core CA RMSD: [measure rmsd $fit2 $fit1]"

  $fit1 delete
  $fit2 delete
  $ox2_all delete

  _common_display
  return [list $ox1 $ox2]
}

proc scene_ox1 {} {
  global OX1_CHAIN OX1_LIG OX1_POS_2_61 OX1_POS_3_33

  set molid 0
  _clear_reps $molid
  _reset_labels
  _common_display

  mol addrep $molid
  mol modselect 0 $molid [_base_receptor_selection $OX1_CHAIN]
  mol modstyle 0 $molid NewCartoon
  mol modcolor 0 $molid ColorID 7
  mol modmaterial 0 $molid AOChalky

  mol addrep $molid
  mol modselect 1 $molid "chain $OX1_CHAIN and resname $OX1_LIG"
  mol modstyle 1 $molid Licorice 0.18 16 16
  mol modcolor 1 $molid ColorID 4
  mol modmaterial 1 $molid AOChalky

  mol addrep $molid
  mol modselect 2 $molid [_water_shell_selection $OX1_CHAIN $OX1_LIG 4.0]
  mol modstyle 2 $molid VDW 0.45 14
  mol modcolor 2 $molid ColorID 0
  mol modmaterial 2 $molid AOChalky

  mol addrep $molid
  mol modselect 3 $molid "protein and chain $OX1_CHAIN and (resid $OX1_POS_2_61 or resid $OX1_POS_3_33)"
  mol modstyle 3 $molid VDW 0.85 16
  mol modcolor 3 $molid ColorID 1
  mol modmaterial 3 $molid AOChalky

  _add_site_labels $molid $OX1_CHAIN [list $OX1_POS_2_61 $OX1_POS_3_33]
  mol showrep $molid 0 on
  display resetview
  scale by 1.15
}

proc scene_ox2 {} {
  global OX2_CHAIN OX2_LIG OX2_POS_2_61 OX2_POS_3_33

  set molid 1
  _clear_reps $molid
  _reset_labels
  _common_display

  mol addrep $molid
  mol modselect 0 $molid [_base_receptor_selection $OX2_CHAIN]
  mol modstyle 0 $molid NewCartoon
  mol modcolor 0 $molid ColorID 3
  mol modmaterial 0 $molid AOChalky

  mol addrep $molid
  mol modselect 1 $molid "chain $OX2_CHAIN and resname $OX2_LIG"
  mol modstyle 1 $molid Licorice 0.18 16 16
  mol modcolor 1 $molid ColorID 1
  mol modmaterial 1 $molid AOChalky

  mol addrep $molid
  mol modselect 2 $molid [_water_shell_selection $OX2_CHAIN $OX2_LIG 4.0]
  mol modstyle 2 $molid VDW 0.45 14
  mol modcolor 2 $molid ColorID 0
  mol modmaterial 2 $molid AOChalky

  mol addrep $molid
  mol modselect 3 $molid "protein and chain $OX2_CHAIN and (resid $OX2_POS_2_61 or resid $OX2_POS_3_33)"
  mol modstyle 3 $molid VDW 0.85 16
  mol modcolor 3 $molid ColorID 11
  mol modmaterial 3 $molid AOChalky

  _add_site_labels $molid $OX2_CHAIN [list $OX2_POS_2_61 $OX2_POS_3_33]
  mol showrep $molid 0 on
  display resetview
  scale by 1.15
}

proc scene_overlay {} {
  global OX1_CHAIN OX2_CHAIN OX1_LIG OX2_LIG
  global OX1_POS_2_61 OX1_POS_3_33 OX2_POS_2_61 OX2_POS_3_33

  _clear_reps 0
  _clear_reps 1
  _reset_labels
  _common_display

  mol addrep 0
  mol modselect 0 0 [_base_receptor_selection $OX1_CHAIN]
  mol modstyle 0 0 NewCartoon
  mol modcolor 0 0 ColorID 7
  mol modmaterial 0 0 AOChalky

  mol addrep 1
  mol modselect 1 0 [_base_receptor_selection $OX1_CHAIN]
  mol modstyle 1 0 Surf 1.4 0.0
  mol modcolor 1 0 ColorID 7
  mol modmaterial 1 0 Transparent

  mol addrep 1
  mol modselect 0 1 [_base_receptor_selection $OX2_CHAIN]
  mol modstyle 0 1 NewCartoon
  mol modcolor 0 1 ColorID 3
  mol modmaterial 0 1 AOChalky

  mol addrep 1
  mol modselect 1 1 [_base_receptor_selection $OX2_CHAIN]
  mol modstyle 1 1 Surf 1.4 0.0
  mol modcolor 1 1 ColorID 3
  mol modmaterial 1 1 Transparent

  mol addrep 0
  mol modselect 2 0 "chain $OX1_CHAIN and resname $OX1_LIG"
  mol modstyle 2 0 Licorice 0.18 16 16
  mol modcolor 2 0 ColorID 4
  mol modmaterial 2 0 AOChalky

  mol addrep 1
  mol modselect 2 1 "chain $OX2_CHAIN and resname $OX2_LIG"
  mol modstyle 2 1 Licorice 0.18 16 16
  mol modcolor 2 1 ColorID 1
  mol modmaterial 2 1 AOChalky

  mol addrep 0
  mol modselect 3 0 [_water_shell_selection $OX1_CHAIN $OX1_LIG 3.5]
  mol modstyle 3 0 VDW 0.42 14
  mol modcolor 3 0 ColorID 0
  mol modmaterial 3 0 AOChalky

  mol addrep 1
  mol modselect 3 1 [_water_shell_selection $OX2_CHAIN $OX2_LIG 3.5]
  mol modstyle 3 1 VDW 0.42 14
  mol modcolor 3 1 ColorID 10
  mol modmaterial 3 1 AOChalky

  mol addrep 0
  mol modselect 4 0 "protein and chain $OX1_CHAIN and (resid $OX1_POS_2_61 or resid $OX1_POS_3_33)"
  mol modstyle 4 0 VDW 0.85 16
  mol modcolor 4 0 ColorID 1
  mol modmaterial 4 0 AOChalky

  mol addrep 1
  mol modselect 4 1 "protein and chain $OX2_CHAIN and (resid $OX2_POS_2_61 or resid $OX2_POS_3_33)"
  mol modstyle 4 1 VDW 0.85 16
  mol modcolor 4 1 ColorID 11
  mol modmaterial 4 1 AOChalky

  _add_site_labels 0 $OX1_CHAIN [list $OX1_POS_2_61 $OX1_POS_3_33]
  _add_site_labels 1 $OX2_CHAIN [list $OX2_POS_2_61 $OX2_POS_3_33]

  display resetview
  scale by 1.10
}

proc render_help {} {
  puts "Suggested workflow:"
  puts "  1. source orexin_vis.tcl"
  puts "  2. orexin_load"
  puts "  3. scene_ox1      ;# render OX1R"
  puts "  4. scene_ox2      ;# render OX2R"
  puts "  5. scene_overlay  ;# render aligned comparison"
  puts ""
  puts "Suggested render settings:"
  puts "  display ambientocclusion on"
  puts "  display shadows on"
  puts "  render TachyonInternal output.tga"
  puts ""
  puts "Presentation angle:"
  puts "  OX1R and OX2R share a highly conserved orthosteric pocket."
  puts "  Selective OX1R antagonism is interesting because OX1R is more tied"
  puts "  to reward/compulsive/addiction-related circuitry, but structure-based"
  puts "  selectivity is hard because the pocket differences are subtle."
}

render_help
