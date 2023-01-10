(TeX-add-style-hook
 "PDE"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz" "crop" "convert={density=300,outext=.png}" "border=0.0cm" "width=18cm" "height=7cm")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "pgfplots"
    "amsmath"
    "physics"))
 :latex)
