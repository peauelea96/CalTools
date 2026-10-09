# Publishing CalTools on GitHub

The project is prepared locally; it has not been published to an external account.

1. Replace `maintainer@example.org` in DESCRIPTION with the real maintainer email,
   complete the author's name and add any coauthors. The address is a placeholder.
2. Review the proposed GPL (>= 2) license and credits in inst/doc/AUDIT.md.
   The supplied prototype declares an adaptation of sampling::calib (GPL >= 2).
3. Run `R CMD build CalTools`, then
   `R CMD check --no-manual CalTools_0.2.0.tar.gz`.
4. Create a repository in your GitHub account and upload the contents of CalTools.
5. Add the actual repository URL and BugReports URL to DESCRIPTION.
6. The included .github/workflows/R-CMD-check.yaml checks future revisions on
   Linux, Windows and macOS. Those remote runs have not yet been executed.

Once published:

```r
# install.packages("remotes")
remotes::install_github("YOUR_USERNAME/CalTools")
```

CalTools is a proposed name; availability on CRAN has not been established.
This is a development release, not a CRAN publication.
