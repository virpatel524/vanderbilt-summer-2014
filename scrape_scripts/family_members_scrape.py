# given a family name, this script pulls all of the family members from mirviewer
__author__ = 'virpatel'
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

    fle = open("mirna_families.txt","r")
    mirna = fle.readlines()
    print mirna
    holder = []
    for i in mirna:
        holder.append(i.split("\n")[0])
    main2homo = {}

    mirna = holder[:]

    for i in mirna[:]:
        main2homo[i] = scrape(i)
        print mirna.index(i)



    return main2homo






def mirnagath():
    mirnas = []
    scraped = scrapeall()
    for key in scraped:
        for allspecies in scraped[key]:
            mirnas.append(allspecies[0][1])
    return mirnas





def exporter():
    mirnas = mirnagath()
    fle = open("homomirna.txt","w")
    for i in mirnas:
        fle.write(i+"\n")

exporter()