# Create the output directory first
#dir.create("docs/walkthrough", recursive = TRUE)

# Check pandoc is available and get its path
rmarkdown::pandoc_exec()

# Use the full pandoc path explicitly
pandoc <- rmarkdown::pandoc_exec()
system(paste0('"', pandoc, '" walkthrough.docx -o docs/walkthrough/index.html --self-contained'))

# Preview html file to edit it
browseURL("docs/walkthrough/index.html")

# If it's too narrow add this code
# Adjust 1000px up or down to change
#May need to remove duplicates of this code
html <- readLines("docs/walkthrough/index.html")

style <- c(
  "<style>",
  "  body { max-width: 1000px; margin: 0 auto; padding: 2em; font-size: 16px; line-height: 1.6; }",
  "</style>"
)

# Insert before closing </head>
head_line <- which(grepl("</head>", html))
html <- append(html, style, after = head_line - 1)

writeLines(html, "docs/walkthrough/index.html")
browseURL("docs/walkthrough/index.html")

#When it looks good, just commit and push edited index.html to github 
