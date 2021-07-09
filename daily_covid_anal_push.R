library(RPushbullet)
tstamp <- format(Sys.time(), "%Y-%m-%d")
pbSetup(apikey = "o.KasCi9GJjdN4Q5301FFeem2lKGdVXfCu",
        conffile = ".rpushbullet.json",
        defdev = "Samsung SM-G977N")

pbPost("note", title = paste0(tstamp, " daily pf completed!"))
