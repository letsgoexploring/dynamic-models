$(document).ready(function () {

    $('#ParameterInput').on('submit', function (e) {
        e.preventDefault();

        // Load parameters
        var periods = parseInt($('#periods').val());
        var alpha = parseFloat($('#alpha').val());
        var delta = parseFloat($('#delta').val());
        var beta = parseFloat($('#beta').val());
        var sigma = parseFloat($('#sigma').val());
        var phi = parseFloat($('#phi').val());
        var rhoa = parseFloat($('#rhoa').val());
        var siga = parseFloat($('#siga').val());

        //  Compute steady state
        var A = 1
        var k_l_ratio = (alpha*A/(Math.pow(beta,-1)+delta-1))**(1/(1-alpha))
        var L = 1 - phi/((1-alpha)*A*Math.pow(k_l_ratio,alpha))
        var K = k_l_ratio*L
        var Y = A*Math.pow(K,alpha)*Math.pow(L,(1-alpha))
        var I = delta*K
        var C = Y - I
        
        // document.write("Steady state:","<br><br>")
        // document.write("A: ",A,"<br>")
        // document.write("K: ",K,"<br>")
        // document.write("C: ",C,"<br>")
        // document.write("L: ",L,"<br>")
        // document.write("I: ",I,"<br>")
        // document.write("Y: ",Y,"<br>")
        // document.write("<br><br>")

        // Construct the coefficients for the log-linearized model
        var phi_1 = -(alpha+(1-alpha)*L)/(1-L)
        var phi_2 = (alpha-1)*(1-beta*(1-delta))
        var phi_3 = (1-alpha)*(1-beta*(1-delta))
        var phi_4 = -rhoa*(1-beta*(1-delta))
        var phi_5 = A*Math.pow(K,alpha-1)*Math.pow(L,1-alpha)
        var phi_6 = alpha*A*Math.pow(K,alpha-1)*Math.pow(L,1-alpha) + 1 - delta
        var phi_7 = -C/K
        var phi_8 = (1-alpha)*A*Math.pow(K,alpha-1)*Math.pow(L,1-alpha)

        // document.write("Coefficients of linearized model:","<br><br>")
        // document.write("phi_1: ",phi_1,"<br>")
        // document.write("phi_2: ",phi_2,"<br>")
        // document.write("phi_3: ",phi_3,"<br>")
        // document.write("phi_4: ",phi_4,"<br>")
        // document.write("phi_5: ",phi_5,"<br>")
        // document.write("phi_6: ",phi_6,"<br>")
        // document.write("phi_7: ",phi_7,"<br>")
        // document.write("phi_8: ",phi_8,"<br>")
        // document.write("<br><br>")
        
        // Compute the solution for the log-linearized model
        var pi_5 = -1/phi_1
        var pi_6 = -alpha/phi_1

        var a = -sigma*phi_7
        var b = phi_2*phi_7  +sigma- sigma*phi_6-sigma*phi_8*pi_6 + phi_3*pi_6*phi_7
        var c = phi_2*phi_6+phi_2*phi_8*pi_6+phi_3*pi_6*(phi_6+phi_8*pi_6)
        var pi_4 = (-b+Math.sqrt(b**2-4*a*c))/2/a

        var top = phi_4 - phi_2*(phi_5+phi_8*pi_5)+sigma*pi_4*(phi_5+phi_8*pi_5) - phi_3*rhoa*pi_5-phi_3*pi_6*(phi_5+phi_8*pi_5)
        var bottom = phi_2*phi_7-sigma*rhoa -sigma*pi_4*phi_7 +phi_3*pi_6*phi_7+sigma
        var pi_3 = top/bottom
        
        var pi_1 = phi_5 + phi_7*pi_3 + phi_8*pi_5
        var pi_2 = phi_6 + phi_7*pi_4 + phi_8*pi_6

        var pi_7 = Y/I*(1 + (1-alpha)*pi_5) - C/I*pi_3
        var pi_8 = Y/I*(alpha + (1-alpha)*pi_6) - C/I*pi_4
        var pi_9 = 1 + (1-alpha)*pi_5
        var pi_10 = alpha + (1-alpha)*pi_6

        // document.write("Solution coefficients:","<br><br>")
        // document.write("pi_1: ",pi_1,"<br>")
        // document.write("pi_2: ",pi_2,"<br>")
        // document.write("pi_3: ",pi_3,"<br>")
        // document.write("pi_4: ",pi_4,"<br>")
        // document.write("pi_5: ",pi_5,"<br>")
        // document.write("pi_6: ",pi_6,"<br>")
        // document.write("pi_7: ",pi_7,"<br>")
        // document.write("pi_8: ",pi_8,"<br>")
        // document.write("pi_9: ",pi_9,"<br>")
        // document.write("pi_10: ",pi_10,"<br>")

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

        var a = 0;
        var k = 0;
        var a_series = [0];
        var k_series = [0];
        var c_series = [0];
        var l_series = [0];
        var i_series = [0];
        var y_series = [0];

        var input = document.getElementById ("stochSim");
        
            if (input.checked == true) {
                for (i = 0; i <= periods; i++) {

                    eps = z.nextGaussian()
                    
                    k = pi_1*a+pi_2*k
                    a = rhoa*a+siga*eps

                    a_series.push(a);
                    k_series.push(k)
                    c_series.push(pi_3*a+pi_4*k);
                    l_series.push(pi_5*a+pi_6*k);
                    i_series.push(pi_7*a+pi_8*k);
                    y_series.push(pi_9*a+pi_10*k);

                }
            } else {
                for (i = 0; i <= periods; i++) {
                    if (i<1) {
                        eps = 1
                    }
                    else {
                        eps=0
                    }

                k = pi_1*a+pi_2*k
                a = rhoa*a+siga*eps

                a_series.push(a);
                k_series.push(k)
                c_series.push(pi_3*a+pi_4*k);
                l_series.push(pi_5*a+pi_6*k);
                i_series.push(pi_7*a+pi_8*k);
                y_series.push(pi_9*a+pi_10*k);

                }
            };

        // document.write("Simulation data:","<br><br>")
        // document.write("a: ",a_series,"<br>")
        // document.write("k: ",k_series,"<br>")
        // document.write("c: ",c_series,"<br>")
        // document.write("l: ",l_series,"<br>")
        // document.write("i: ",i_series,"<br>")
        // document.write("y: ",y_series,"<br>")

        sessionStorage.setItem('a_series',JSON.stringify(a_series));
        sessionStorage.setItem('k_series',JSON.stringify(k_series));
        sessionStorage.setItem('c_series',JSON.stringify(c_series));
        sessionStorage.setItem('y_series',JSON.stringify(y_series));
        sessionStorage.setItem('l_series',JSON.stringify(l_series));
        sessionStorage.setItem('i_series',JSON.stringify(i_series));
        sessionStorage.setItem('periods',periods);


        $('#aProcess').highcharts({
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
                text: 'TFP',
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
                minRange: 0.0001,
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
                name: 'a',
                // color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*a_series[i])/100000,
                            // y: g_series[i],
                        })
                    }
                    return data;
                })()
            }
            ]
        });

        $('#kProcess').highcharts({
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
                text: 'Capital',
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
                minRange: 0.0001,
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
                name: 'k',
                // color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*k_series[i])/100000,
                        })
                    }
                    return data;
                })()
            }
            ]
        });

        $('#cProcess').highcharts({
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
                text: 'Consumption',
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
                minRange: 0.0001,
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
                name: 'c',
                // color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*c_series[i])/100000,
                        })
                    }
                    return data;
                })()
            }
            ]
        });
        $('#yProcess').highcharts({
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
                text: 'Output',
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
                minRange: 0.0001,
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
                name: 'y',
                // color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*y_series[i])/100000,
                        })
                    }
                    return data;
                })()
            }
            ]
        });
        $('#lProcess').highcharts({
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
                text: 'Labor',
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
                minRange: 0.0001,
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
                name: 'l',
                // color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*l_series[i])/100000,
                        })
                    }
                    return data;
                })()
            }
            ]
        });
        $('#iProcess').highcharts({
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
                text: 'Investment',
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
                minRange: 0.0001,
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
                name: 'i',
                // color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*i_series[i])/100000,
                        })
                    }
                    return data;
                })()
            }
            ]
        });
    }) //.trigger('submit');
});

function reloadFunction() {
    sessionStorage.clear();
    location.reload();
}

function downloadFunction() {
    
    if (sessionStorage.getItem('a_series') == null){
        window.alert('Run the simulation first.')
    } else {
        var a_data = JSON.parse(sessionStorage.getItem('a_series'));
        var k_data = JSON.parse(sessionStorage.getItem('k_series'));
        var c_data = JSON.parse(sessionStorage.getItem('c_series'));
        var l_data = JSON.parse(sessionStorage.getItem('l_series'));
        var y_data = JSON.parse(sessionStorage.getItem('y_series'));
        var i_data = JSON.parse(sessionStorage.getItem('i_series'));
        var periods = parseInt(sessionStorage.getItem('periods'));

        var row1 = [];
        row1.push("Period");
        row1.push("A");
        row1.push("K");
        row1.push("C");
        row1.push("Y");
        row1.push("L");
        row1.push("I");

        var rows = [row1]
        for (i = 0; i <= periods; i++) {
            rows.push([i,a_data[i],k_data[i],c_data[i],y_data[i],l_data[i],i_data[i]]);
        }

        let csvContent = "data:text/csv;charset=utf-8,";
        rows.forEach(function(rowArray){
           let row = rowArray.join(",");
           csvContent += row + "\r\n"; // add carriage return
        }); 

        var encodedUri = encodeURI(csvContent);
        var link = document.createElement("a");
        link.setAttribute("href", encodedUri);
        link.setAttribute("download", "rbc_data.csv");
        document.body.appendChild(link); // Required for FF

        link.click(); // This will download the data file named "rbc_data.csv".
    }

}