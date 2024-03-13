library("readxl")

scores <- read_excel("tests/scores.xlsx")

countries <- read_csv("tests/iso-en.csv")
name_to_iso <- countries$iso_alpha_3
names(name_to_iso) <- countries$name
name_to_iso["Netherlands"] <- "NLD"
name_to_iso["Slovak Republic"] <- "SVK"
name_to_iso["Türkiye"] <- "TUR"
name_to_iso["Korea"] <- "KOR"

scores$country <- str_replace(scores$country, "\\*", "")
scores$country_iso <- name_to_iso[scores$country]

write_csv(scores, "tests/scores.csv")
