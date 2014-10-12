# same scraper as before but this time it uses the threshold definitions to determine

import urllib
import re
from scrapy.selector import Selector

def scrape(name):

    bigfile = []
    eachspecies  = []







    urltest = "http://people.csail.mit.edu/akiezun/microRNAviewer/all_" + name + ".html"

    page = urllib.urlopen(urltest)
    page = page.read()



    sel = Selector(text=page)
    holder1 = []
    holder2 = []

    table = sel.xpath("//tr")
    table =  table[4:]
    selectors = []

    for i in table:
        selectors.append(Selector(text=i.extract()))


    for i in selectors:
        holder1.append(i.xpath("//tr/td[@title]").extract())


    rawdata = holder1[:]


    allspecies = []
    for i in rawdata:
        minilst = []
        for var in i:

            var =  str(var[var.find("title=\"")+7:var.find("\"",var.find("title")+8)])



            var =  re.sub("\([^>]+\)", " ", var)
            var = var.split(" ")

            temp1 = var[0]
            temp2 = var[-1].split(":")[-1]

            temp3 = temp1[:3]
            box = [str(temp3),str(temp1),str(temp2)]
            minilst.append(box)
        allspecies.append(minilst)

    return allspecies




def scrapeall():

    fle = open("mirna.txt","r")
    mirna = fle.readlines()
    holder = []
    for i in mirna:
        holder.append(i.split("\n")[0])
    main2homo = {}

    mirna = holder[:]

    for i in mirna[:]:
        main2homo[i] = scrape(i)
        print mirna.index(i)


    return main2homo




scrapeddata = scrapeall()

def threshold(scrapedata,num):
    homfam = []
    for key in scrapeddata:

        for homo in scrapeddata[key]:
            famem = ""
            for species in homo:
                if float(species[-1]) >= float(num) and species[0]!="xtr":
                    if species!= homo[-1]:
                        famem += species[0]+"|miRviewer:"+species[1]+" "
                    else:
                        famem += species[0]+"|miRviewer:"+species[1]
            homfam.append(famem)
    return homfam



def filewrite(scrapeddata):
    for i in [.5,.6,.75,.8,.85,.95,1.0]:
        temp = threshold(scrapeddata,i)
        fle = open("fams/"+str(int(i*100))+"fams.txt","w")

        for i in temp:
            if i!= "":
                fle.write(i+"\n")



filewrite(scrapeddata)

