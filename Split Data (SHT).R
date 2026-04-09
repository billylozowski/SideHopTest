# split.txt separates the master .txt into its individual trials
# since each participant already has an assigned code, extract this and the trial from the file name
# - i.e. "MBB001_THEIA" + :SHL_001" becomes "MBB001_THEIA_SHL_001.csv"

library(tidyverse)
library(readr)

# define the function
split.txt <- function(path,
                      out.dir = "Data/Split Data",
                      delim = "\t",
                      connection.size = 10485760) {  # 10 MB
   
   # read in the master data file
   Sys.setenv("VROOM_CONNECTION_SIZE" = connection.size)
   
   d <- read_delim(path, col_names = FALSE, 
                   delim = "\t", # files are separated by tabs
                   show_col_types = FALSE, na = character(), 
                   progress = FALSE)
   
   # extract the first row which contains the full file name
   trial.names <- d %>%
      slice(1) %>%
      select(-1) %>%
      unlist(use.names = FALSE) %>%
      as.character() %>%
      trimws()
   
   # replace missing or blank trial names
   if (any(is.na(trial.names) | trial.names == "")) {
      blanks <- which(is.na(trial.names) | trial.names == "")
      trial.names[blanks] <- paste0("trial_", blanks)
   }
   
   # group column indices by trial name
   idx.by.trial <- split(seq_along(trial.names), trial.names)
   
   # create a safe file name for each trial
   safe.name <- function(x) {
      x %>%
         str_replace_all("[^-_.A-Za-z0-9]+", "-") %>% # make filesystem-safe
         str_replace("^[-_]+|[-_]+$", "")             # trim leading/trailing - or _
   }
   
   # create the output directory if it doesn't already exist
   if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
   # ensure the "Data" folder exists, then save the list of trial names
   if (!dir.exists("Data")) dir.create("Data", recursive = TRUE)
   
   # extract participant folder and trial file from full path
   extract_trial_parts <- function(x) {
      parts <- str_split(x, "[/\\\\]+", simplify = TRUE)                # split by "/" or "\"
      parts <- parts[parts != ""]                                       # remove any blank elements
      n <- length(parts)                                                # count the parts
      tibble(participant = if (n >= 2) parts[n - 1] else NA_character_, # second-to-last element (participant code)
             trial.file  = if (n >= 1) parts[n]     else NA_character_) # last element (file name)
   }
   
   # collector for the output summary (list of trial names/files)
   trial.summary <- tibble(original.file = character(),
                           participant   = character(),
                           trial.file    = character(),
                           trial.stem    = character(),
                           trial.name    = character(),
                           file.path     = character())
   
   # loop through each trial
   iwalk(
      .x = idx.by.trial,
      .f = \(idx, nm) {
         # always keep the first column "ITEM" (this is the sample number)
         # idx positions are relative to columns AFTER dropping the first;
         # add +1 to map them back to absolute positions in 'd'
         cols <- c(1L, idx + 1L)
         
         # drop the first row (trial name labels)
         out.df <- d[-1, cols, drop = FALSE]
         
         # split full path into participant folder and trial file
         parts       <- extract_trial_parts(nm)
         participant <- parts$participant %>%
            str_remove("_THEIA$")
         trial.file  <- parts$trial.file
         
         # remove original extension only (i.e. ".c3d")
         trial.stem <- str_remove(trial.file, "\\.[^.]*$")
         
         # final trial file name: "<participant>_<trial stem>"
         trial.name <- str_c(participant, "_", trial.stem)
         out.path   <- file.path(out.dir, paste0(safe.name(trial.name), ".csv"))
         
         # save the individual trial .csv
         write_csv(out.df, out.path, col_names = FALSE)
         
         # append to the summary table
         trial.summary <<- bind_rows(
            trial.summary,
            tibble(original.file = nm,
                   participant   = participant,
                   trial.file    = trial.file,
                   trial.stem    = trial.stem,
                   trial.name    = trial.name,
                   file.path     = out.path)
         )
      }
   )
   
   # save the complete trial summary df
   write_csv(trial.summary %>%
                distinct(trial.name, .keep_all = TRUE) %>%
                arrange(participant, trial.stem),
             file.path("Data", "Trial Names.csv"))
   
   # return the summary invisibly for program use
   invisible(trial.summary)
}

# call the function
split.txt("Data/Exported SHT Data/MBB 2025_09.txt")
split.txt("Data/Exported SHT Data/WBB 2025_09.txt")
split.txt("Data/Exported SHT Data/WVB 2025_08.txt")
