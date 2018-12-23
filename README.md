# Plubell_2017_PAW
### Data from Plubell et al., 2017 processed with PAW pipeline

A lot of work has happened since the 2017 MCP paper. In that first publication, IRS was done in Excel. It really is a simple enough idea to do that. It is labor intensive and possibly error-prone, though. I have a [data analysis pipeline](https://github.com/pwilmart/PAW_pipeline.git) originally developed for SEQUEST search results that is written in Python. It has nice protein inference and protein grouping steps. It has been updated to work with [Comet](http://comet-ms.sourceforge.net/) search results (at least up through 2016 versions).

One thing I do not like about Proteome Discoverer, is the protein inference and how shared peptides are used in quantification. The PSM export files from PD have all of the same fields that [the PAW pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2816815/) needs, along with the reporter ion information. It is possible to take the confidently identified PSMs in the PD exports and run them through the later stages of the PAW pipeline. Support for TMT reporter ions was added.

The PAW pipeline used MSConvert of the [ProteoWizard](http://proteowizard.sourceforge.net/) package to extract the MS2 scan information for Comet searches. Support to extract the reporter ion scan peak heights was added. Now the PAW pipeline can take TMT data exported from PD and produce protein-level quantitative reports, or data straight from RAW files in a full open source pipeline.

The data from the original publication was depositied in [PRIDE](https://www.ebi.ac.uk/pride/archive/projects/PXD005953) and has been re-analyzed with Proteome Discoverer 2.2, PAW/Comet, and MaxQuant. Notebooks for analysis of these different workflows will be added eventually.

One very important part of doing an IRS experiment that uses pooled internal standards, is making sure that those channels are correctly specified. There is nothing about the IRS procedure that has any knowledge of the correct channels to use for the internal standards except you! If you make a mistake in the standard channel designations, your data will get messed up. Like most computer use, there is no real way to protect you from yourself. The quality and accuracy of the record keeping is on you.

That said, we can actually get the computers to help us double check our records of which channels were the pooled standard channels. The "auto_finder_PAW" notebook show you how to see which channels are the most similar in a TMT plex without specifying any sample information. The notebook reads the PAW results files, but the concepts would apply to other results files (PD or MaxQuant).

I will add more content to this repository as time allows.

December 23, 2018 - Phil W.
