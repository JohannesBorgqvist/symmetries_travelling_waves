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
    "./Input/RD"
    "./Input/TV_l_1"
    "article"
    "art12"
    "geometry"
    "parskip"
    "inputenc"
    "amsmath"
    "amssymb"
    "amsfonts"
    "amsthm")
   (LaTeX-add-amsthm-newtheorems
    "theorem"))
 :latex)

