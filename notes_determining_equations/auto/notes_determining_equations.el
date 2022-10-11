(TeX-add-style-hook
 "notes_determining_equations"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper" "12pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=0.7in") ("parskip" "parfill") ("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "./Input/intro"
    "./Input/det_eq"
    "./Input/det_eq_1"
    "./Input/det_eq_2"
    "./Input/det_eq_3"
    "./Input/det_eq_3_l_0"
    "./Input/det_eq_3_l_1"
    "./Input/det_eq_3_l_2"
    "article"
    "art12"
    "geometry"
    "parskip"
    "inputenc"
    "amsmath"
    "amssymb"
    "amsfonts"
    "amsthm"))
 :latex)

