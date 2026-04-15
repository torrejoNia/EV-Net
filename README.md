
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- github markdown built using
rmarkdown::render("README.Rmd",output_format = "md_document")
-->

<p align="center">

<img src="man/figures/logo.png" width="200"/>
</p>

**EV-Net: A computational framework to model extracellular
vesicles-mediated communication.** EV-Net enables the identification and
prioritization of EV cargo molecules with high regulatory potential in a
recipient tissue of interest.

## Installation

EV-Net is an R package that requires R version 4.3.2 or higher for
installation.

Installation time depends on the number of dependencies already
installed. In our testing, it ranged from approximately 5 to 20 minutes.
You can install EV-Net (and its required dependencies) directly from
GitHub:

``` r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}

devtools::install_github("torrejoNia/EV-Net")

library(EVNet)
```

EV-Net has been tested on Windows, Linux and MacOS. Most recently tested
with R version 4.5.2.

## Learn more about EV-Net

For detailed documentation, tutorials, and examples, visit the: 👉
[EV-Net page](https://torrejonia.github.io/EV-Net/)

## Authors

[Estefania
Torrejón](https://scholar.google.com/citations?user=c5RvnlUAAAAJ&hl=en&oi=ao),
[Joeri
Sleegers](https://scholar.google.com/citations?user=hRSv5kwAAAAJ&hl=en&oi=ao),
[Rune
Matthiesen](https://scholar.google.com/citations?user=tsjn1BcAAAAJ&hl=en&oi=ao),
[Maria Paula
Macedo](https://scholar.google.com/citations?user=qjZpKj0AAAAJ&hl=en&oi=ao),
[Anaïs
Baudot](https://scholar.google.com/citations?user=rzoYM1cAAAAJ&hl=en&oi=ao),
[Rita Machado de
Oliveira](https://scholar.google.com/citations?user=b55XP-oAAAAJ&hl=en&oi=ao)

## References

**Pre-print:** Torrejon, E., Sleegers, J., Matthiesen, R., Macedo, M.
P., Baudot, A., & Machado de Oliveira, R. (2026). EV-Net: A
computational framework to model extracellular vesicles-mediated
communication. bioRxiv, 2026-04.

Read it
[here](https://www.biorxiv.org/content/10.64898/2026.04.02.716053v1)

## Funding

The core of this work was developed during the EMBO scientific exchange
grant 10987 carried out by Estefania Torrejón in Anaïs Baudot lab.
Additional work was supported by the Research Unit iNOVA4Health
(UID/4462/2025) and by the Associated Laboratory LS4FUTURE
(LA/P/0087/2020), both financially supported by Fundação para a Ciência
e Tecnologia / Ministério da Educação, Ciência e Inovação, and the EVCA
Twinning Project (Horizon GA n° 101079264) financially supported by the
European Union.

The aims of this study contribute to the ERDERA project (participant:
Anaïs Baudot), which has received funding from the European Union’s
Horizon Europe research and innovation programme under grant agreement
N°101156595. Anaïs Baudot was also supported by France 2030 state
funding managed by the National Research Agency with the reference
“ANR-22-PESN-0013”. Rita Machado de Oliveira was supported by Fundação
para a Ciência e a Tecnologia (2022.05764.PTDC). Maria Paula Macedo was
supported by the European Commission CORDIS Pas Gras Project
(101080329). Estefania Torrejón’s PhD fellowship is funded by the
iNOVA4Health-FCT fellowship UI/BD/154345/2022.

Presentations in conferences were supported by the European Union
Twining Project EVCA (Horizon GA n° 101079264) and International Society
of Computational Biology Student Council (ISCB-SC).

### If you have problems installing dependencies

Some systems may fail to install `ComplexHeatmap` and `Seurat`
automatically when installing EV-Net, because they require additional
repositories or dependencies.

Please install them manually first:

``` r
# Install BiocManager if needed
install.packages("BiocManager") 
# Install ComplexHeatmap from Bioconductor
BiocManager::install("ComplexHeatmap") 

# Install Seurat and SeuratObject from CRAN
install.packages(c("SeuratObject", "Seurat"))
```

After that, install **EV-Net**:

``` r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}

devtools::install_github("torrejoNia/EV-Net")
```
