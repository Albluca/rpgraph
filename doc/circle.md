This document describes how to produce a circular principal graph.

Building a principal circle
---------------------------

The `simple_circle` dataset included in the package describe points placed on a three dimensional circle and We can use it to test the usage of principal cicrles

``` r
library(rpgraph)
```

    ## Loading required package: rJava

    ## 
    ## Attaching package: 'rpgraph'

    ## The following object is masked from 'package:base':
    ## 
    ##     Filter

``` r
Data <- simple_circle
Results <- computeElasticPrincipalGraph(Data = Data, NumNodes = 40, Method = 'CircleConfiguration')
```

    ## Configuring engine ......[1] "Empty initialization"
    ## [1] ""
    ## [1] "Running engine"

Now `Results` will be a list and the first element will contain the processed principal graph. Diagnostic information can be obtained using

``` r
plotMSDEnergyPlot(Results[[1]], Main = "Pincipal Circle", Cex.Main = 1)
```

![](circle_files/figure-markdown_github/unnamed-chunk-2-1.png)

and

``` r
accuracyComplexityPlot(Results[[1]], Main = "Pincipal Circle", Cex.Main = 1, Mode = 5)
```

![](circle_files/figure-markdown_github/unnamed-chunk-3-1.png)

It it also possible to zoom into a specific area of the accuracy/complexity plot by using the Xlims parameter.

``` r
accuracyComplexityPlot(Results[[1]], Main = "Pincipal Circle", Cex.Main = 1, Xlims = c(.97, .98))
```

![](circle_files/figure-markdown_github/unnamed-chunk-4-1.png)

Data can be plotted in 2D using the R built-in functions

``` r
plotData2D(Data = simple_circle, PrintGraph = Results[[1]],
           GroupsLab = rep(1, nrow(simple_circle)), Xlab = "Dimension 1", Ylab = "Dimension 2")
```

![](circle_files/figure-markdown_github/unnamed-chunk-5-1.png)

or plotly, which produces an interactive plot. Using plotly interactivelly requires running the code in RStudio (does it?)

``` r
plotData2D(Data = simple_circle, PrintGraph = Results[[1]], Plot.ly = TRUE,
           GroupsLab = factor(rep(1, nrow(simple_circle))), Xlab = "Dimension 1", Ylab = "Dimension 2")
```

    ## Warning in arrange_impl(.data, dots): '.Random.seed' n'est pas un vecteur
    ## d'entiers, mais est de type 'NULL', et sera donc ignoré

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

    ## Warning: plotly.js doesn't yet support line.width arrays, track this issue for progress
    ##   https://github.com/plotly/plotly.js/issues/147

<!--html_preserve-->

<script type="application/json" data-for="htmlwidget-87178b981fa7ec816a57">{"x":{"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"xaxis":{"domain":[0,1],"title":"Dimension 1"},"yaxis":{"domain":[0,1],"title":"Dimension 2"},"title":"","hovermode":"closest"},"source":"A","config":{"modeBarButtonsToAdd":[{"name":"Collaborate","icon":{"width":1000,"ascent":500,"descent":-50,"path":"M487 375c7-10 9-23 5-36l-79-259c-3-12-11-23-22-31-11-8-22-12-35-12l-263 0c-15 0-29 5-43 15-13 10-23 23-28 37-5 13-5 25-1 37 0 0 0 3 1 7 1 5 1 8 1 11 0 2 0 4-1 6 0 3-1 5-1 6 1 2 2 4 3 6 1 2 2 4 4 6 2 3 4 5 5 7 5 7 9 16 13 26 4 10 7 19 9 26 0 2 0 5 0 9-1 4-1 6 0 8 0 2 2 5 4 8 3 3 5 5 5 7 4 6 8 15 12 26 4 11 7 19 7 26 1 1 0 4 0 9-1 4-1 7 0 8 1 2 3 5 6 8 4 4 6 6 6 7 4 5 8 13 13 24 4 11 7 20 7 28 1 1 0 4 0 7-1 3-1 6-1 7 0 2 1 4 3 6 1 1 3 4 5 6 2 3 3 5 5 6 1 2 3 5 4 9 2 3 3 7 5 10 1 3 2 6 4 10 2 4 4 7 6 9 2 3 4 5 7 7 3 2 7 3 11 3 3 0 8 0 13-1l0-1c7 2 12 2 14 2l218 0c14 0 25-5 32-16 8-10 10-23 6-37l-79-259c-7-22-13-37-20-43-7-7-19-10-37-10l-248 0c-5 0-9-2-11-5-2-3-2-7 0-12 4-13 18-20 41-20l264 0c5 0 10 2 16 5 5 3 8 6 10 11l85 282c2 5 2 10 2 17 7-3 13-7 17-13z m-304 0c-1-3-1-5 0-7 1-1 3-2 6-2l174 0c2 0 4 1 7 2 2 2 4 4 5 7l6 18c0 3 0 5-1 7-1 1-3 2-6 2l-173 0c-3 0-5-1-8-2-2-2-4-4-4-7z m-24-73c-1-3-1-5 0-7 2-2 3-2 6-2l174 0c2 0 5 0 7 2 3 2 4 4 5 7l6 18c1 2 0 5-1 6-1 2-3 3-5 3l-174 0c-3 0-5-1-7-3-3-1-4-4-5-6z"},"click":"function(gd) { \n        // is this being viewed in RStudio?\n        if (location.search == '?viewer_pane=1') {\n          alert('To learn about plotly for collaboration, visit:\\n https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html');\n        } else {\n          window.open('https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html', '_blank');\n        }\n      }"}],"modeBarButtonsToRemove":["sendDataToCloud"]},"data":[{"x":[0.344771792,-0.955250694,1.125562866,1.123355382,-0.764858395,0.468984683,0.166718467,1.11605487,-0.423529056,-0.681326311,1.217838221,0.261524331,0.630954897,-0.55212558,1.037008662,0.132124062,-0.189259728,-0.635592209,-0.375107762,-0.046345795,-0.951346202,1.050970803,1.105837912,0.14601591,1.14428442,-0.440476345,0.471888495,0.108087053,0.934099918,0.722244296,0.833009658,0.81987548,1.257197588,0.563777867,-0.004311245,0.659723309,-0.655837104,-0.06831261,0.86868235,-0.610435358,-0.750596642,0.391028274,-0.334636613,1.48667494,1.057282403,0.989070115,0.122659032,-0.719179898,-0.021925146,-0.961460601,-0.215846202,0.815771373,-0.523199588,-0.561903242,-0.354606903,-0.934533555,0.590323062,-0.8376083,-0.437409227,0.018442466,-0.330124523,0.904439432,1.295207993,-0.190061282,-0.411494546,0.548869381,-0.364383161,0.002652376,0.791717468,0.468494723,0.697996006,1.576519202,0.548588305,0.430703015,0.995790922,-0.616429529,0.974379886,0.333020671,-0.396384235,-0.556621786,0.368640195,0.395514809,1.059069255,0.475287148,0.202031865,1.017409475,0.885511771,1.169915209,-0.059835898,-0.497417472,0.79634477,0.945219104,-0.338336158,0.948814039,0.511490152,-0.521917348,0.327658673,1.512521924,0.439262928,0.739264127,-0.840855703,-0.39134672,-0.049046337,-0.311076097,-0.784809328,1.436018891,-0.297790811,-0.844839958,0.823578411,0.855972251,1.383642784,0.377126398,-0.513642739,0.682884004,-0.543298193,-0.336399284,0.286914609,-0.470023038,0.643280381,0.644248941,0.23204318,1.144273826,-0.760663064,-0.669773804,-0.046753843,0.646198734,-0.80020262,0.834039854,-0.168252494,-0.357227753,0.039616814,-0.195588311,0.627048974,0.924182333,0.639545391,-0.932420144,-0.518491535,-0.09398354,-0.926609859,0.661248226,-0.385469027,-0.953927182,-0.337944381,1.064784006,0.982691345,-0.850237613,0.485233835,-0.560720625,1.147182685,-0.545865868,0.268374573,0.919695044,0.604768817,0.490524679,1.07165315,-0.740625168,0.234567834,1.097480588,0.500303003,-0.681746833,-0.243037609,0.296122123,-0.171687269,-0.931571509,0.933828826,0.034872447,-0.277226583,-0.574423757,1.505073107,-0.772232971,0.178329473,-0.720912222,-0.072572097,0.158025228,1.369560785,0.551062142,0.921465295,0.284836035,0.674985162,0.545459372,0.884265875,0.788731123,0.160924689,0.42482261,1.096412266,-0.539337595,-0.812660824,1.067849121,1.142761261,1.133440624,1.210884581,0.641257805,-0.26952982,0.89155228,-0.498951207,0.632129847,-0.771755135,-0.573582406,1.295583423,0.537931085],"y":[1.251424192,0.172530751,1.171294192,0.118630795,0.20375007,1.196843876,-0.91946646,0.338881482,-0.690766056,0.369034533,-0.313711909,-0.85616991,0.886782093,-0.543105909,0.748697735,1.039766741,-0.519926014,0.225542127,1.438904471,1.553928266,0.155989014,0.069425255,0.073259959,1.333294693,-0.055213236,0.93316242,1.254289875,-0.731413391,1.420319929,1.016087095,-0.531676095,-0.592007673,0.678743336,1.010337145,1.179456206,1.292704666,-0.431389277,-0.032123594,-0.365869976,-0.524308372,0.803626265,-0.502698917,1.135869303,0.540499086,-0.427675172,-0.189523395,-0.907596463,-0.079263209,1.2729387,0.028641347,0.807985118,-0.428484402,-0.510301902,-0.39742945,-0.349402465,0.444834156,-0.783374436,0.07483333,-0.538697038,-0.49634986,1.008363234,-0.17412716,-0.285003033,1.228744656,1.44053122,-0.671639387,0.786305098,-0.832677078,1.343399451,1.097322881,-0.52042826,0.402661255,1.078182645,-0.866832595,0.511359545,-0.567700592,0.515614226,1.62261841,-0.152940394,0.948943086,1.297343334,-0.713565366,1.312239866,-0.991036133,1.040883331,0.856852544,0.863334106,0.248819346,-0.757704524,0.943512883,-0.326121024,0.551888095,-0.530125651,0.479302439,-0.674516117,0.46837848,1.34611853,0.375546179,1.134995154,-0.791824265,0.564940123,-0.60275349,-0.593658409,-0.641600716,-0.052383772,0.330791265,1.023286601,0.02165745,-0.68722582,-0.554226168,0.463016579,-0.826059225,-0.193802059,1.334151995,-0.547328561,1.391743756,-0.954509138,1.080139229,-0.572756195,0.785819692,-0.656686487,0.712950105,-0.37051404,0.795079238,1.607422611,-0.553567338,0.346883827,1.065776555,-0.900476668,-0.874632883,-0.979909299,0.543650123,-0.594759645,-0.129866366,1.075025502,0.327119436,-0.604333847,1.208447252,-0.13284487,1.36466079,1.313021938,0.584316188,-0.705658956,0.330681778,0.585068768,-0.106474705,0.991538367,-0.383882592,-0.051331406,-0.754732787,1.289380512,-0.520750873,1.06851635,-0.674259053,0.66860621,0.199374535,1.278230387,-0.625616149,-0.673850784,0.06336231,0.056042261,-0.524443543,-0.275051224,0.410239474,0.743818333,1.043410401,-0.489275272,0.639627439,0.140110841,0.276632738,-0.318865814,0.794855703,-0.511222468,-0.7564861,0.017927075,1.373205378,1.241413726,1.402892285,0.879294701,1.357406877,-0.407844247,1.063265893,-0.607600778,-0.955604118,1.225515344,-0.405202784,-0.231138964,0.001111294,0.662253695,-0.101729537,-0.074034287,-0.731837811,-0.793209719,-0.495554555,0.207642151,-0.66985694,-0.386261847,0.649517683,0.036290237,1.101616643],"mode":"markers","text":["R_1","R_2","R_3","R_4","R_5","R_6","R_7","R_8","R_9","R_10","R_11","R_12","R_13","R_14","R_15","R_16","R_17","R_18","R_19","R_20","R_21","R_22","R_23","R_24","R_25","R_26","R_27","R_28","R_29","R_30","R_31","R_32","R_33","R_34","R_35","R_36","R_37","R_38","R_39","R_40","R_41","R_42","R_43","R_44","R_45","R_46","R_47","R_48","R_49","R_50","R_51","R_52","R_53","R_54","R_55","R_56","R_57","R_58","R_59","R_60","R_61","R_62","R_63","R_64","R_65","R_66","R_67","R_68","R_69","R_70","R_71","R_72","R_73","R_74","R_75","R_76","R_77","R_78","R_79","R_80","R_81","R_82","R_83","R_84","R_85","R_86","R_87","R_88","R_89","R_90","R_91","R_92","R_93","R_94","R_95","R_96","R_97","R_98","R_99","R_100","R_101","R_102","R_103","R_104","R_105","R_106","R_107","R_108","R_109","R_110","R_111","R_112","R_113","R_114","R_115","R_116","R_117","R_118","R_119","R_120","R_121","R_122","R_123","R_124","R_125","R_126","R_127","R_128","R_129","R_130","R_131","R_132","R_133","R_134","R_135","R_136","R_137","R_138","R_139","R_140","R_141","R_142","R_143","R_144","R_145","R_146","R_147","R_148","R_149","R_150","R_151","R_152","R_153","R_154","R_155","R_156","R_157","R_158","R_159","R_160","R_161","R_162","R_163","R_164","R_165","R_166","R_167","R_168","R_169","R_170","R_171","R_172","R_173","R_174","R_175","R_176","R_177","R_178","R_179","R_180","R_181","R_182","R_183","R_184","R_185","R_186","R_187","R_188","R_189","R_190","R_191","R_192","R_193","R_194","R_195","R_196","R_197","R_198","R_199","R_200"],"hoverinfo":"text","type":"scatter","name":"1","marker":{"size":10,"sizemode":"area","fillcolor":"rgba(255,0,0,0.5)","color":"rgba(255,0,0,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y"},{"x":[0.810866832733154,1.26050508022308,-0.0782072991132736,-0.369355380535126,-0.527670443058014,-0.363566726446152,0.0919722467660904,1.25129771232605,-0.00503998063504696,0.329349339008331,1.07398772239685,1.05758595466614,-0.615729510784149,-0.787837088108063,0.666174948215485,-0.726640164852142,0.783772885799408,0.657215058803558,-0.431185454130173,-0.625630140304565,0.525259613990784,0.961667239665985,-0.724139750003815,1.15993094444275,-0.28298944234848,0.500220656394958,0.394703805446625,-0.782622933387756,1.30288624763489,-0.242885559797287,0.161249592900276,0.874017000198364,-0.795489728450775,0.989025712013245,1.16949808597565,0.249100863933563,-0.459748268127441,-0.509968876838684,-0.155808508396149,0.893154799938202],"y":[1.09959101676941,0.400779455900192,-0.698876202106476,-0.588553845882416,-0.353705108165741,1.09206175804138,-0.739714324474335,0.110751658678055,1.29317426681519,1.27986836433411,-0.196634978055954,0.631725013256073,-0.224024266004562,0.374673366546631,1.18651127815247,0.532939314842224,-0.559764981269836,-0.643291473388672,0.941240251064301,0.678119659423828,-0.710587561130524,0.777801632881165,-0.0934766605496407,-0.0408112332224846,1.22401940822601,1.23466169834137,-0.770296633243561,0.0532983466982841,0.265615105628967,-0.660125553607941,1.28276491165161,0.937544405460358,0.211207941174507,-0.336552023887634,0.519335925579071,-0.771389663219452,-0.483004122972488,0.800857484340668,1.27614772319794,-0.455843508243561],"mode":"markers","text":["V_1","V_2","V_3","V_4","V_5","V_6","V_7","V_8","V_9","V_10","V_11","V_12","V_13","V_14","V_15","V_16","V_17","V_18","V_19","V_20","V_21","V_22","V_23","V_24","V_25","V_26","V_27","V_28","V_29","V_30","V_31","V_32","V_33","V_34","V_35","V_36","V_37","V_38","V_39","V_40"],"hoverinfo":"text","type":"scatter","name":"Graph","marker":{"size":10,"sizemode":"area","fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)","line":{"color":"transparent"}},"xaxis":"x","yaxis":"y"},{"x":[0.0919722467660904,-0.0782072991132736],"y":[-0.739714324474335,-0.698876202106476],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.527670443058014,-0.615729510784149],"y":[-0.353705108165741,-0.224024266004562],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.666174948215485,0.810866832733154],"y":[1.18651127815247,1.09959101676941],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.787837088108063,-0.726640164852142],"y":[0.374673366546631,0.532939314842224],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.783772885799408,0.657215058803558],"y":[-0.559764981269836,-0.643291473388672],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.431185454130173,-0.363566726446152],"y":[0.941240251064301,1.09206175804138],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.726640164852142,-0.625630140304565],"y":[0.532939314842224,0.678119659423828],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.657215058803558,0.525259613990784],"y":[-0.643291473388672,-0.710587561130524],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.961667239665985,1.05758595466614],"y":[0.777801632881165,0.631725013256073],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.615729510784149,-0.724139750003815],"y":[-0.224024266004562,-0.0934766605496407],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[1.25129771232605,1.15993094444275],"y":[0.110751658678055,-0.0408112332224846],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[1.15993094444275,1.07398772239685],"y":[-0.0408112332224846,-0.196634978055954],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.363566726446152,-0.28298944234848],"y":[1.09206175804138,1.22401940822601],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.329349339008331,0.500220656394958],"y":[1.27986836433411,1.23466169834137],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.500220656394958,0.666174948215485],"y":[1.23466169834137,1.18651127815247],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.525259613990784,0.394703805446625],"y":[-0.710587561130524,-0.770296633243561],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.724139750003815,-0.782622933387756],"y":[-0.0934766605496407,0.0532983466982841],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[1.26050508022308,1.30288624763489],"y":[0.400779455900192,0.265615105628967],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[1.30288624763489,1.25129771232605],"y":[0.265615105628967,0.110751658678055],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.0782072991132736,-0.242885559797287],"y":[-0.698876202106476,-0.660125553607941],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.242885559797287,-0.369355380535126],"y":[-0.660125553607941,-0.588553845882416],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.00503998063504696,0.161249592900276],"y":[1.29317426681519,1.28276491165161],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.161249592900276,0.329349339008331],"y":[1.28276491165161,1.27986836433411],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.810866832733154,0.874017000198364],"y":[1.09959101676941,0.937544405460358],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.874017000198364,0.961667239665985],"y":[0.937544405460358,0.777801632881165],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.782622933387756,-0.795489728450775],"y":[0.0532983466982841,0.211207941174507],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.795489728450775,-0.787837088108063],"y":[0.211207941174507,0.374673366546631],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[1.07398772239685,0.989025712013245],"y":[-0.196634978055954,-0.336552023887634],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[1.05758595466614,1.16949808597565],"y":[0.631725013256073,0.519335925579071],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[1.16949808597565,1.26050508022308],"y":[0.519335925579071,0.400779455900192],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.394703805446625,0.249100863933563],"y":[-0.770296633243561,-0.771389663219452],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.249100863933563,0.0919722467660904],"y":[-0.771389663219452,-0.739714324474335],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.369355380535126,-0.459748268127441],"y":[-0.588553845882416,-0.483004122972488],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.459748268127441,-0.527670443058014],"y":[-0.483004122972488,-0.353705108165741],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.625630140304565,-0.509968876838684],"y":[0.678119659423828,0.800857484340668],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.509968876838684,-0.431185454130173],"y":[0.800857484340668,0.941240251064301],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.28298944234848,-0.155808508396149],"y":[1.22401940822601,1.27614772319794],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[-0.155808508396149,-0.00503998063504696],"y":[1.27614772319794,1.29317426681519],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.989025712013245,0.893154799938202],"y":[-0.336552023887634,-0.455843508243561],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"},{"x":[0.893154799938202,0.783772885799408],"y":[-0.455843508243561,-0.559764981269836],"mode":"lines","text":["",""],"hoverinfo":"text","type":"scatter","showlegend":false,"name":"Graph","line":{"width":1,"fillcolor":"rgba(0,0,0,0.5)","color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y"}],"base_url":"https://plot.ly"},"evals":["config.modeBarButtonsToAdd.0.click"],"jsHooks":[]}</script>
<!--/html_preserve-->
This commands will produce a list of warnings, which can be ignored. Unfortunately there is not an easy way to fix it at this time.

The plotly graph can be exported on the web, for example on [plot.ly](http://plot.ly) using the instruction provided [here](http://plot.ly/r/getting-started/).

Data can also be plotted in 3D using the functionalities provided by the `rgl` package

``` r
plotData3D(Data = simple_circle, PrintGraph = Results[[1]], Plot.ly = FALSE,
           GroupsLab = factor(rep(1, nrow(simple_circle))), NodeSizeMult = 0.05,
           Xlab = "Dimension 1", Ylab = "Dimension 2", Zlab = "Dimension 3")
```

For `rgl` to work correctly on MacOS, a proper X11 environment need to be installed. The most common solution is the installation of [XQuartz](http://www.xquartz.org/). If the code crashes R, try using

``` r
library("rgl")
open3d()
```

before invoking the `plotData3D` function and/or rebooting the system.

It is also possible to produce 3D plots using plot.ly.

``` r
plotData3D(Data = simple_circle, PrintGraph = Results[[1]], Plot.ly = TRUE,
           GroupsLab = factor(rep(1, nrow(simple_circle))),
           Xlab = "Dimension 1", Ylab = "Dimension 2", Zlab = "Dimension 3")
```

Similarly to 2D, it is possible to export the plot to a web resource.

To see how different populations distribute among the nodes of the graph by providing a population identifier to each point in the original data (The `Categories` vector). In the following example, three populations are randomly assigned.

``` r
Net <- ConstructGraph(Results = Results[[1]], DirectionMat = NULL, Thr = 0.05)
TaxonList <- getTaxonMap(Results = Results[[1]], Data = Data)

InfoData <- plotPieNet(Results = Results[[1]], Data = simple_circle, NodeSizeMult = 4,
                       Categories = factor(sample(1:3, nrow(simple_circle), replace = TRUE)),
           Graph = Net, TaxonList = TaxonList, LayOut = 'circle', Main = "Pincipal Circle")
```
