(TeX-add-style-hook
 "individual_figures_initial_conditions"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz" "crop" "convert={density=200,outext=.png}" "border=0.4cm")))
   (TeX-run-style-hooks
    "latex2e"
    "./Input/u"
    "standalone"
    "standalone10"
    "pgfplots"
    "amsmath"
    "physics"
    "xcolor")
   (LaTeX-add-xcolor-definecolors
    "exp_1"
    "exp_2"
    "log_1"
    "log_2"
    "gen_1"
    "gen_2"))
 :latex)

