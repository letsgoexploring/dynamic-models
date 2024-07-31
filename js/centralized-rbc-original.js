// Functions

// Define a rounding function for display precision
round = function(x,n) { return Math.round(Math.pow(10,n)*x)/Math.pow(10,n); };

// Define a square function for root-finding
sqr = function(x) { return x*x; };

// Do the charting things
function build_charts(allData) {
    var periods = allData.t.length;
    $('#eProc').highcharts({
        chart: {
            type: 'line',
            marginRight: 10,
        },

        credits: {
            enabled: false
        },

        plotOptions: {
            series: {
                lineWidth: 3,
                marker : {
                    enabled : false,
                    radius : 3},
                animation: {
                    duration: 10000     //Controls the time required for plot to be fully rendered.
                }
            }
        },

        title: {
            text: 'TFP shock process',
            useHTML:true,
            style: {
                "fontSize": "15px"
            }
        },
        xAxis: {
            text: 'Time Period'
        },
        yAxis: {
            title: {
                text: ''
            },
            labels: {
                formatter: function()
                {
                    return round(this.value,2);
                }
            },
            minRange: 0,
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
            csv: {
                columnHeaderFormatter: function (series, key) {
                    if (key === 'y') {
                        return series.name;
                    } else {
                        return 'Time';
                    }
                }
            },
        },
        series: [{
            name: 'epsilon',
            data: (function () {
                var data = [];
                for (i = 0; i <= periods; i++) {
                    data.push({
                        x: allData.t[i],
                        y: round(allData.e[i],5),
                    })
                }

                return data;
            })()
        },
        ]
    });

    $('#aProc').highcharts({
        chart: {
            type: 'line',
            marginRight: 10,
        },

        credits: {
            enabled: false
        },

        plotOptions: {
            series: {
                lineWidth: 3,
                marker : {
                    enabled : false,
                    radius : 3},
                animation: {
                    duration: 10000     //Controls the time required for plot to be fully rendered.
                }
            }
        },

        title: {
            text: 'TFP',
            useHTML:true,
            style: {
                "fontSize": "15px"
            }
        },
        xAxis: {
            text: 'Time Period'
        },
        yAxis: {
            title: {
                text: ''
            },
            labels: {
                formatter: function()
                {
                    return round(this.value,2);
                }
            },
            minRange: 0,
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
            csv: {
                columnHeaderFormatter: function (series, key) {
                    if (key === 'y') {
                        return series.name;
                    } else {
                        return 'Time';
                    }
                }
            },
        },
        series: [{
            name: 'a',
            data: (function () {
                var data = [];
                for (i = 0; i <= periods; i++) {
                    data.push({
                        x: allData.t[i],
                        y: round(allData.a[i],5),
                    })
                }

                return data;
            })()
        },
        ]
    });

    $('#lProc').highcharts({
        chart: {
            type: 'line',
            marginRight: 10,
        },

        credits: {
            enabled: false
        },

        plotOptions: {
            series: {
                lineWidth: 3,
                marker : {
                    enabled : false,
                    radius : 3},
                animation: {
                    duration: 10000     //Controls the time required for plot to be fully rendered.
                }
            }
        },

        title: {
            text: 'Labor',
            useHTML:true,
            style: {
                "fontSize": "15px"
            }
        },
        xAxis: {
            text: 'Time Period'
        },
        yAxis: {
            title: {
                text: ''
            },
            labels: {
                formatter: function()
                {
                    return round(this.value,2);
                }
            },
            minRange: 0,
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
            csv: {
                columnHeaderFormatter: function (series, key) {
                    if (key === 'y') {
                        return series.name;
                    } else {
                        return 'Time';
                    }
                }
            },
        },
        series: [{
            name: 'l',
            data: (function () {
                var data = [];
                for (i = 0; i <= periods; i++) {
                    data.push({
                        x: allData.t[i],
                        y: round(allData.l[i],5),
                    })
                }

                return data;
            })()
        },
        ]
    });

    $('#yProc').highcharts({
        chart: {
            type: 'line',
            marginRight: 10,
        },

        credits: {
            enabled: false
        },

        plotOptions: {
            series: {
                lineWidth: 3,
                marker : {
                    enabled : false,
                    radius : 3},
                animation: {
                    duration: 10000     //Controls the time required for plot to be fully rendered.
                }
            }
        },

        title: {
            text: 'Output',
            useHTML:true,
            style: {
                "fontSize": "15px"
            }
        },
        xAxis: {
            text: 'Time Period'
        },
        yAxis: {
            title: {
                text: ''
            },
            labels: {
                formatter: function()
                {
                    return round(this.value,2);
                }
            },
            minRange: 0,
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
            csv: {
                columnHeaderFormatter: function (series, key) {
                    if (key === 'y') {
                        return series.name;
                    } else {
                        return 'Time';
                    }
                }
            },
        },
        series: [{
            name: 'y',
            data: (function () {
                var data = [];
                for (i = 0; i <= periods; i++) {
                    data.push({
                        x: allData.t[i],
                        y: round(allData.y[i],5),
                    })
                }

                return data;
            })()
        },
        ]
    });

    $('#kProc').highcharts({
        chart: {
            type: 'line',
            marginRight: 10,
        },

        credits: {
            enabled: false
        },

        plotOptions: {
            series: {
                lineWidth: 3,
                marker : {
                    enabled : false,
                    radius : 3},
                animation: {
                    duration: 10000     //Controls the time required for plot to be fully rendered.
                }
            }
        },

        title: {
            text: 'Capital',
            useHTML:true,
            style: {
                "fontSize": "15px"
            }
        },
        xAxis: {
            text: 'Time Period'
        },
        yAxis: {
            title: {
                text: ''
            },
            labels: {
                formatter: function()
                {
                    return round(this.value,2);
                }
            },
            minRange: 0,
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
            csv: {
                columnHeaderFormatter: function (series, key) {
                    if (key === 'y') {
                        return series.name;
                    } else {
                        return 'Time';
                    }
                }
            },
        },
        series: [{
            name: 'k',
            data: (function () {
                var data = [];
                for (i = 0; i <= periods; i++) {
                    data.push({
                        x: allData.t[i],
                        y: round(allData.k[i],5),
                    })
                }

                return data;
            })()
        },
        ]
    });

    $('#iProc').highcharts({
        chart: {
            type: 'line',
            marginRight: 10,
        },

        credits: {
            enabled: false
        },

        plotOptions: {
            series: {
                lineWidth: 3,
                marker : {
                    enabled : false,
                    radius : 3},
                animation: {
                    duration: 10000     //Controls the time required for plot to be fully rendered.
                }
            }
        },

        title: {
            text: 'Investment',
            useHTML:true,
            style: {
                "fontSize": "15px"
            }
        },
        xAxis: {
            text: 'Time Period'
        },
        yAxis: {
            title: {
                text: ''
            },
            labels: {
                formatter: function()
                {
                    return round(this.value,2);
                }
            },
            minRange: 0,
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
            csv: {
                columnHeaderFormatter: function (series, key) {
                    if (key === 'y') {
                        return series.name;
                    } else {
                        return 'Time';
                    }
                }
            },
        },
        series: [{
            name: 'Investment',
            data: (function () {
                var data = [];
                for (i = 0; i <= periods; i++) {
                    data.push({
                        x: allData.t[i],
                        y: round(allData.i[i],5),
                    })
                }

                return data;
            })()
        },
        ]
    });

    $('#cProc').highcharts({
        chart: {
            type: 'line',
            marginRight: 10,
        },

        credits: {
            enabled: false
        },

        plotOptions: {
            series: {
                lineWidth: 3,
                marker : {
                    enabled : false,
                    radius : 3},
                animation: {
                    duration: 10000     //Controls the time required for plot to be fully rendered.
                }
            }
        },

        title: {
            text: 'Consumption',
            useHTML:true,
            style: {
                "fontSize": "15px"
            }
        },
        xAxis: {
            text: 'Time Period'
        },
        yAxis: {
            title: {
                text: ''
            },
            labels: {
                formatter: function()
                {
                    return round(this.value,2);
                }
            },
            minRange: 0,
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
            csv: {
                columnHeaderFormatter: function (series, key) {
                    if (key === 'y') {
                        return series.name;
                    } else {
                        return 'Time';
                    }
                }
            },
        },
        series: [{
            name: 'Consumption',
            data: (function () {
                var data = [];
                for (i = 0; i <= periods; i++) {
                    data.push({
                        x: allData.t[i],
                        y: round(allData.c[i],5),
                    })
                }

                return data;
            })()
        },
        ]
    });
}


$(document).ready(function () {

    $('#ParameterInput').on('submit', function (e) {
        e.preventDefault();
        var reqData = {
            A: parseFloat($('#tfp').val()),
            alpha: parseFloat($('#alpha').val()),
            delta: parseFloat($('#delta').val()),
            beta: parseFloat($('#beta').val()),
            sigma: parseFloat($('#sigma').val()),
            eta: parseFloat($('#eta').val()),
            phi: parseFloat($('#phi').val()),
            rhoa: parseFloat($('#rhoA').val()),
            sige: parseFloat($('#sigA').val()),
            periods: parseInt($('#periods').val()),
            stochSim: $('input[name=simType]:checked').val() == 'stoch' ? 1 : 0
        }

        $.ajax(
            'https://letsgoexploring.techb.us/api/v1/centralized-rbc-with-labor-simulation/',
            {
                dataType: 'json',
                data: reqData
            }
        ).done(function(data) {
            build_charts(data);
            sessionStorage.setItem('t',JSON.stringify(data.t));
            sessionStorage.setItem('e',JSON.stringify(data.e));
            sessionStorage.setItem('a',JSON.stringify(data.a));
            sessionStorage.setItem('k',JSON.stringify(data.k));
            sessionStorage.setItem('c',JSON.stringify(data.c));
            sessionStorage.setItem('i',JSON.stringify(data.i));
            sessionStorage.setItem('y',JSON.stringify(data.y));
            sessionStorage.setItem('l',JSON.stringify(data.l));
            
        }).error(function() {
            alert('Something has gone wrong. Please try again.')
        });
    })
});

function reloadFunction() {
    sessionStorage.clear();
    location.reload();
}

function downloadFunction() {
    
    if (sessionStorage.getItem('k') == null){
        window.alert('Run the simulation first.')
    } else {
        var t = JSON.parse(sessionStorage.getItem('t'));
        var periods = t[t.length-1]
        var e = JSON.parse(sessionStorage.getItem('e'));
        var a = JSON.parse(sessionStorage.getItem('a'));
        var k = JSON.parse(sessionStorage.getItem('k'));
        var c = JSON.parse(sessionStorage.getItem('c'));
        var invest = JSON.parse(sessionStorage.getItem('i'));
        var y = JSON.parse(sessionStorage.getItem('y'));
        var l = JSON.parse(sessionStorage.getItem('l'));

        var row1 = [];
        row1.push("period");
        row1.push("productivity shock");
        row1.push("tfp");
        row1.push("capital");
        row1.push("output");
        row1.push("consumption");
        row1.push("investment");
        row1.push("labor");


        var rows = [row1];
        for (i = 0; i <= periods; i++) {
            rows.push([i,e[i],a[i],k[i],y[i],c[i],invest[i],l[i]])
        }

        let csvContent = "data:text/csv;charset=utf-8,";
        rows.forEach(function(rowArray){
           let row = rowArray.join(",");
           csvContent += row + "\r\n"; // add carriage return
        }); 

        var encodedUri = encodeURI(csvContent);
        var link = document.createElement("a");
        link.setAttribute("href", encodedUri);
        link.setAttribute("download", "centralized-rbc-with-labor-data.csv");
        document.body.appendChild(link); // Required for FF

        link.click(); // This will download the data file named "centralized-rbc-with-labor-data.csv".
    }

}
