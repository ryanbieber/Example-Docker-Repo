# Docker Build Model Template Repo

## Step 1:
### Make sure you have docker installed on your desktop

## Step 2:
### Place your functions in functions.R, your hardcoded values you want in parameters.R(Not sensitive ones like DB credentials), place your actually modelling script inside model.R.

## Step 3: 
### Use the setup.R file to install the packages you are going to use, I have the timeseries package I built as an example package.

## Step 4:
### Make sure to save everything then make your shell towards the folder with dockerfile in it. After you are there, 
### you should run the command `docker build -t model-image .`

### This will build the image locally and it may take a while since it has to download a lot of stuff and it ends up  being around 700 mB. After this is done running you then want to run `docker run -it model-image:latest` This will run the image from the entrypoint you have laid out which is the main.R script. If you are looking at your shell I used the -it tag to denote interactivity. This is nice if you want to see if work and it should by default when the image is done running delete itself.


## Step 5:
### To deploy this image onto any cloud deployment solution is fairly trivial as you are just going to place this image in a cluster of your choice and they are already many guides out there for that. I just wanted to make a simple example of how much simplier it is to deploy things in containers than onto your local machine. You can schedule this model to run using a cronjob, make a REST API, the possibilities are endless. Have fun turning everything you own into containers!
 
