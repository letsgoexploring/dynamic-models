$(document).ready(function () {


    $('#ParameterInput').on('submit', function (e) {
        e.preventDefault();
        var periods = parseInt($('#periods').val());
        var muX = parseFloat($('#muX').val());
        var rhoX = parseFloat($('#rhoX').val());
        var sigEps = parseFloat($('#sigEps').val());
        var x0 = parseFloat($('#x0').val());

        var interval = 0;
        var enableMarks = true;


        var z = new ziggurat();


        if (periods<=30) {
            interval = 5;
        }
        else if (30<periods&&periods<=100) {
            interval = 20;
            enableMarks = false;
        }
        else if (100<periods&&periods<=500) {
            interval = 100;
            enableMarks = false;
        }
        else {
            interval = 500;
            enableMarks = false;
        };



        var xSeries = [x0];
        var eSeries = [0];

        var input = document.getElementById ("stochSim");
            if (input.checked == true) {
                for (i = 0; i <= periods; i++) {
                    e = z.nextGaussian()
                    xSeries.push(rhoX*xSeries[xSeries.length-1]+sigEps*e);
                    eSeries.push(e);

                }
            } else {
                for (i = 0; i <= periods; i++) {
                    if (i<1) {
                        xSeries.push(rhoX*xSeries[xSeries.length-1]+sigEps);
                        eSeries.push(sigEps);
                    }
                    else {
                        xSeries.push(rhoX*xSeries[xSeries.length-1]);
                        eSeries.push(0);
                    }

                }
            };

        sessionStorage.setItem('eSeries',JSON.stringify(eSeries));
        sessionStorage.setItem('xSeries',JSON.stringify(xSeries));
        sessionStorage.setItem('periods',periods);


        $('#shockProcess').highcharts({
            chart: {
                type: 'line',
                marginRight: 10,
            },

            credits: {
                enabled: false
            },

            plotOptions: {
                series: {
                    lineWidth: 2,
                    marker : {
                        enabled : false,
                        radius : 3},
                    animation: {
                        duration: 10000     //Controls the time required for plot to be fully rendered.
                    }
                }
            },

            title: {
                text: 'Shock Process',
                style: {
                    "fontSize": "15px" 
                }
            },
            xAxis: {
                // type: 'datetime',
                // tickPixelInterval: 150,
                tickInterval: interval,
                text: 'Time'
            },
            yAxis: {
                title: {
                    text: ''
                },
                labels: {
                    formatter: function() 
                    {
                        return Math.round(this.value*100)/100;
                    }
                },
                minRange: 1,
                plotLines: [{
                    value: 0,
                    width: 2,
                    color: '#808080'
                }]
            },

            legend: {
                enabled: false
            },
            exporting: {
                enabled: true,
                buttons: {
                    contextButton: {
                        menuItems: ['downloadPNG']
                    },
                },
            },
            series: [{
                name: 'e',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*eSeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            ]
        });

        $('#ar1Process').highcharts({
            chart: {
                type: 'line',
                marginRight: 10,
            },

            credits: {
                enabled: false
            },

            plotOptions: {
                series: {
                    lineWidth: 2,
                    marker : {
                        enabled : false,
                        radius : 3},
                    animation: {
                        duration: 10000     //Controls the time required for plot to be fully rendered.
                    }
                }
            },

            title: {
                text: 'AR(1) Process',
                style: {
                    "fontSize": "15px" 
                }
            },
            xAxis: {
                // type: 'datetime',
                // tickPixelInterval: 150,
                tickInterval: interval,
                text: 'Time'
            },
            yAxis: {
                title: {
                    text: ''
                },
                labels: {
                    formatter: function() 
                    {
                        return Math.round(this.value*100)/100;
                    }
                },
                minRange: 1,
                plotLines: [{
                    value: 0,
                    width: 2,
                    color: '#808080'
                }]
            },

            legend: {
                enabled: false
            },
            exporting: {
                enabled: true,
                buttons: {
                    contextButton: {
                        menuItems: ['downloadPNG']
                    },
                },
            },
            series: [{
                name: 'x',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*xSeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            ]
        });
    }) //.trigger('submit');
});

function reloadFunction() {
    sessionStorage.clear();
    location.reload();
}

function downloadFunction() {
    
    if (sessionStorage.getItem('eSeries') == null){
        window.alert('Run the simulation first.')
    } else {
        var e = JSON.parse(sessionStorage.getItem('eSeries'));
        var x = JSON.parse(sessionStorage.getItem('xSeries'));
        var periods = parseInt(sessionStorage.getItem('periods'));

        var row1 = [];
        row1.push("Period");
        row1.push("Shock process");
        row1.push("AR(1) process");

        var rows = [row1]
        for (i = 0; i <= periods; i++) {
            rows.push([i,e[i],x[i]]);
        }

        let csvContent = "data:text/csv;charset=utf-8,";
        rows.forEach(function(rowArray){
           let row = rowArray.join(",");
           csvContent += row + "\r\n"; // add carriage return
        }); 

        var encodedUri = encodeURI(csvContent);
        var link = document.createElement("a");
        link.setAttribute("href", encodedUri);
        link.setAttribute("download", "ar1_data.csv");
        document.body.appendChild(link); // Required for FF

        link.click(); // This will download the data file named "ar1_data.csv".
    }

}