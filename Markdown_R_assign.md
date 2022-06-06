---
title: "R_assign1.md"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
First I will load all the required packages. 
```{r}
library(janitor)
library(tidyverse)
library(plyr)
library(dplyr )
library(ggplot2)
library(reshape2)
```
After loading the required packages, read the required files into R. 
```{r}
fang<-read_tsv("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt")
view(fang)
snp<-read_tsv("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/snp_position.txt")
view(snp)
```
Now I have both snp and fang files loaded and stored in objects in r.
```{r}
fang
snp
```

This tells me that fang is a 2782X986 tibble and snp is a 983x15 tibble. Considerably smaller
```{r}
object.size(fang) #23.14 MB, or 23124584 bytes
object.size(snp) #359.384 KB, or 359384 bytes
```
So to sum it up, fang is 2782 rows, 986 columns, and 23.14 MB, while snp is 983 rows, 15 columns, and 359.384 KB. 
```{r}
fangcols<-colnames(fang)
snpcols<- colnames(snp)
```
For both files, SNP_ID is the first column. 
```{r}
maize <- filter(fang, `Group` %in% c('ZMMLR','ZMMMR','ZMMIL'))
teosinte <- filter(fang, `Group` %in% c('ZMPBA','ZMPIL','ZMPJA'))
```
This has created subsets of the maize and teosinte values. 
```{r}
view(teosinte)
```

Just checking, it looks good so far. Next I need to transcribe my teosinte and maize files so that they rows and columns run in the other direction. 
```{r}
transteo <- t(teosinte)
transmaize <- t(maize)
```
Next I will cut my SNP file so that I have the columns SNP_ID, Chromosome, and Position. 
```{r}
snpsnip <- select(snp, c('SNP_ID','Chromosome','Position'))
view(snpsnip)
```
Looks good! I need to remove a few unnecessary rows from my files, make headers, and check my dimensions. 
```{r}
trimteo <- row_to_names(transteo, 3, remove_row = TRUE, remove_rows_above = TRUE)
trimmaize <- row_to_names(transmaize, 3, remove_row = TRUE, remove_rows_above = TRUE)
#summary(trimteo) if you want to see the teosinte summary
#summary(snpsnip) if you want to see the trimmed snp file summary
```
My dataframes are equal length and all cleaned up, so now I will join by the common column.
```{r}
teosnp <- cbind(snpsnip,trimteo)
teosintesnp1<-subset(teosnp, Chromosome!="unknown" & Chromosome!="multiple")
teosintesnp<-subset(teosintesnp1, Position!="multiple")
maisnp <- cbind(snpsnip,trimmaize)
maizesnp1<-subset(maisnp, Chromosome!="unknown" & Chromosome!="multiple")
maizesnp<-subset(maizesnp1, Position!="multiple")

```
Now the files are bound together, with the first three columns coming from my snp file.
I have also removed the entries that had unknown and multiple chromosomes. 
I am now ready to move to the next step.
Now I need a file with missing data encoded by "?" and one with missing data encoded by "-".
My file already has missing data encoded by "?". 
Next step is to sort by position values for my existing files. After this point, I will make another file with decreasing position values and "-" in place of missing values. 
```{r}
maizesnp$Position = as.numeric(as.character(maizesnp$Position))
is.numeric(maizesnp$Position)
maize_quest <- maizesnp[order(maizesnp$Position) ,]
teosintesnp$Position = as.numeric(as.character(teosintesnp$Position))
is.numeric(teosintesnp$Position)
teo_quest <- teosintesnp[order(teosintesnp$Position) ,]
```
I had to convert the Position column from character to numeric in order to sort numerically. 
The next step is to reverse sort by Position and replace ? with -.

```{r}
maize_dash1 <- data.frame(lapply(maize_quest, gsub, pattern = "[?]", replacement = "-"))
teo_dash1 <- data.frame(lapply(teo_quest, gsub, pattern = "[?]", replacement = "-"))
maize_dash1$Position = as.numeric(as.character(maize_dash1$Position))
teo_dash1$Position = as.numeric(as.character(teo_dash1$Position))
maize_dash <-maize_dash1[order(rev(maize_dash1$Position)),]
teo_dash <-teo_dash1[order(rev(teo_dash1$Position)),]
maize_quest$Position = as.numeric(as.character(maize_quest$Position))
teo_quest$Position = as.numeric(as.character(teo_quest$Position))
```
Now I'm ready to split into separate files for maize and teosinte data. 
I need chr 1-10 of maize and teo, ascending position (maize_quest, teo_quest) and also chr 1-10 of maize and teo, descending position (maize_dash, teo_dash).
```{r}
maize_increase_1<-subset(maize_quest, Chromosome=="1")
write.csv(maize_increase_1,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_1.csv")
#
maize_increase_2<-subset(maize_quest, Chromosome=="2")
write.csv(maize_increase_2,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_2.csv")
#
maize_increase_3<-subset(maize_quest, Chromosome=="3")
write.csv(maize_increase_3,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_3.csv")
#
maize_increase_4<-subset(maize_quest, Chromosome=="4")
write.csv(maize_increase_4,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_4.csv")
#
maize_increase_5<-subset(maize_quest, Chromosome=="5")
write.csv(maize_increase_5,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_5.csv")
#
maize_increase_6<-subset(maize_quest, Chromosome=="6")
write.csv(maize_increase_6,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_6.csv")
#
maize_increase_7<-subset(maize_quest, Chromosome=="7")
write.csv(maize_increase_7,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_7.csv")
#
maize_increase_8<-subset(maize_quest, Chromosome=="8")
write.csv(maize_increase_8,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_8.csv")
#
maize_increase_9<-subset(maize_quest, Chromosome=="9")
write.csv(maize_increase_9,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_9.csv")
#
maize_increase_10<-subset(maize_quest, Chromosome=="10")
write.csv(maize_increase_10,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_inc_10.csv")
#
teo_increase_1<-subset(teo_quest, Chromosome=="1")
write.csv(teo_increase_1,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_1.csv")
#
teo_increase_2<-subset(teo_quest, Chromosome=="2")
write.csv(teo_increase_2,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_2.csv")
#
teo_increase_3<-subset(teo_quest, Chromosome=="3")
write.csv(teo_increase_3,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_3.csv")
#
teo_increase_4<-subset(teo_quest, Chromosome=="4")
write.csv(teo_increase_4,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_4.csv")
#
teo_increase_5<-subset(teo_quest, Chromosome=="5")
write.csv(teo_increase_5,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_5.csv")
#
teo_increase_6<-subset(teo_quest, Chromosome=="6")
write.csv(teo_increase_6,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_6.csv")
#
teo_increase_7<-subset(teo_quest, Chromosome=="7")
write.csv(teo_increase_7,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_7.csv")
#
teo_increase_8<-subset(teo_quest, Chromosome=="8")
write.csv(teo_increase_8,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_8.csv")
#
teo_increase_9<-subset(teo_quest, Chromosome=="9")
write.csv(teo_increase_9,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_9.csv")
#
teo_increase_10<-subset(teo_quest, Chromosome=="10")
write.csv(teo_increase_10,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_inc_10.csv")
#
maize_decrease_1<-subset(maize_dash, Chromosome=="1")
write.csv(maize_decrease_1,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_1.csv")
#
maize_decrease_2<-subset(maize_dash, Chromosome=="2")
write.csv(maize_decrease_2,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_2.csv")
#
maize_decrease_3<-subset(maize_dash, Chromosome=="3")
write.csv(maize_decrease_3,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_3.csv")
#
maize_decrease_4<-subset(maize_dash, Chromosome=="4")
write.csv(maize_decrease_4,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_4.csv")
#
maize_decrease_5<-subset(maize_dash, Chromosome=="5")
write.csv(maize_decrease_5,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_5.csv")
#
maize_decrease_6<-subset(maize_dash, Chromosome=="6")
write.csv(maize_decrease_6,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_6.csv")
#
maize_decrease_7<-subset(maize_dash, Chromosome=="7")
write.csv(maize_decrease_7,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_7.csv")
#
maize_decrease_8<-subset(maize_dash, Chromosome=="8")
write.csv(maize_decrease_8,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_8.csv")
#
maize_decrease_9<-subset(maize_dash, Chromosome=="9")
write.csv(maize_decrease_9,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_9.csv")
#
maize_decrease_10<-subset(maize_dash, Chromosome=="10")
write.csv(maize_decrease_10,"/Users/gracecarey/Documents/GitHub/r_assign_546//maize_dec_10.csv")
#
teo_decrease_1<-subset(teo_dash, Chromosome=="1")
write.csv(teo_decrease_1,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_1.csv")
#
teo_decrease_2<-subset(teo_dash, Chromosome=="2")
write.csv(teo_decrease_2,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_2.csv")
#
teo_decrease_3<-subset(teo_dash, Chromosome=="3")
write.csv(teo_decrease_3,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_3.csv")
#
teo_decrease_4<-subset(teo_dash, Chromosome=="4")
write.csv(teo_decrease_4,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_4.csv")
#
teo_decrease_5<-subset(teo_dash, Chromosome=="5")
write.csv(teo_decrease_5,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_5.csv")
#
teo_decrease_6<-subset(teo_dash, Chromosome=="6")
write.csv(teo_decrease_6,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_6.csv")
#
teo_decrease_7<-subset(teo_dash, Chromosome=="7")
write.csv(teo_decrease_7,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_7.csv")
#
teo_decrease_8<-subset(teo_dash, Chromosome=="8")
write.csv(teo_decrease_8,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_8.csv")
#
teo_decrease_9<-subset(teo_dash, Chromosome=="9")
write.csv(teo_decrease_9,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_9.csv")
#
teo_decrease_10<-subset(teo_dash, Chromosome=="10")
write.csv(teo_decrease_10,"/Users/gracecarey/Documents/GitHub/r_assign_546//teo_dec_10.csv")
#
```
Now on to visualization!

To do these steps, I will use the "maize_dash" and "teo_dash" files since order does not matter. I choose the dash files because the question mark appears to mess up a lot of R codes.
```{r}
maize_df <- maize_dash
teo_df <- teo_dash
maize_df$Position = as.character(as.double(maize_df$Position))
teo_df$Position = as.character(as.double(teo_df$Position))
maize_df %>% pivot_longer(!Chromosome, names_to = "Position", values_to= "BP", )%>%  {.} -> maize_long
teo_df %>% pivot_longer(!Chromosome, names_to = "Position", values_to= "BP", )%>%  {.} -> teo_long
maize_long$Chromosome = as.numeric(as.character(maize_long$Chromosome))
teo_long$Chromosome = as.numeric(as.character(teo_long$Chromosome))
dim(maize_long)
mutate(maize_long, maize_long$Species <- ("Maize"))
maize_long1 <- mutate(maize_long, maize_long$Species==("Maize"))
maize_long%>% 
  select(Chromosome, Position) %>%
  mutate(
    Species = "Maize"
  )%>% {.} -> maize_mut
view(maize_mut)
teo_long%>% 
  select(Chromosome, Position) %>%
  mutate(
    Species = "Teosinte"
     )%>% {.} -> teo_mut

bind <- bind_rows(maize_mut, teo_mut)
view(bind)
bind%>% 
  select(Chromosome, Position,Species) %>%
  mutate(
    Chromnumber = Chromosome
    )%>% {.} -> bind1
```

Now I have the key file to make visualization1 from!

```{r}
is.numeric(bind1$Chromosome)
bind1$Chromosome = as.character(as.numeric(bind1$Chromosome))
```

And here is the plot which shows distribution of Chromosomes by species: 
```{r}
ggplot(bind1, aes(x=Chromosome, fill= Species, color= Species)) + geom_bar(bins=10, position = "dodge")

```
Now I will move on to step two of the Visualization. I will go back to the steps in my code where I split into maize and teosinte, and rerun through using a combined maize+teo Group selection. I will also bind to the trimmed snp file again. 
```{r}

full_mt <- filter(fang, `Group` %in% c('ZMMLR','ZMMMR','ZMMIL','ZMPBA','ZMPIL','ZMPJA'))
#view(full_mt)
trans_full_mt <- t(full_mt)
trimfull <- row_to_names(trans_full_mt, 3, remove_row = TRUE, remove_rows_above = TRUE)
fullsnp <- cbind(snpsnip,trimfull)
trimfull1<-subset(fullsnp, Chromosome!="unknown" & Chromosome!="multiple")
fullsnpcombo<-subset(trimfull1, Position!="multiple")
#view(fullsnpcombo)
```
Now I have one huge file from which to find homozygous and heterozygous samples. 
```{r}
longest_pivot <- pivot_longer(fullsnpcombo, c(2:2524), names_to = "Sample", values_to = "NT")
longest_pivot_trim <- longest_pivot[-c(1,2),]
head(longest_pivot_trim)
long_df <- longest_pivot_trim
head(long_df)
long_df$Zygosity = ifelse(long_df$NT %in% c("A/A","T/T","C/C","G/G","-/-") ,"Homozygous",
                  ifelse(long_df$NT %in% c("A/A","T/T","C/C","G/G","-/-"), "Heterozygous",
                         "Heterozygous"))                  
Zygosity <- long_df
Zygosity$Group <- str_extract(Zygosity$Sample,"(\\w+)") 
head(Zygosity)
```
At this point, I have labeled cells with heterozygous or homozygous, as appropriate. I also created another column which contains Group IDs. 
I am now ready to visualize this part of the data. 
```{r}
Zy_try <- Zygosity
head(Zy_try)                  
Zy_try$Group <- str_extract(Zy_try$Sample,"(\\w+)") 
head(Zy_try)
Zy_1 <- filter(Zy_try, Group != "Chromosome")
Zy_2 <- filter(Zy_1, Group != "Position")
ggplot(Zy_2, aes(x=Group, fill= Zygosity, color= Zygosity)) + geom_bar(bins=10, position = "dodge")

```

In the above graph, we can see Zygosity by group.

```{r}
Zy_fun <- Zy_2
Zy_fun$Species = ifelse(Zy_fun$Group %in% c('ZMMLR','ZMMMR','ZMMIL') ,"Maize",
                          ifelse(long_df$NT %in% c('ZMPBA','ZMPIL','ZMPJA'), "Teosinte",
                                 "Teosinte"))   
ggplot(Zy_fun, aes(x=Species, fill= Zygosity, color= Zygosity)) + geom_bar(bins=10, position = "dodge")

```
And finally, above see the Zygosity of maize and teosinte, compared. This last graph is my voluntary visualization. 



