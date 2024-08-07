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
        var eta = parseFloat($('#eta').val());
        var rhoa = parseFloat($('#rhoa').val());
        var siga = parseFloat($('#siga').val());

        //  Compute steady state
        var A = 1
        var k_l_ratio = (alpha*A/(Math.pow(beta,-1)+delta-1))**(1/(1-alpha))

        // Solve nonlinear expression for L. References:
        // https://web.archive.org/web/20180821215841/http://numericjs.com/documentation.html
        // https://web.archive.org/web/20180821220006/http://numericjs.com/workshop.php?link=4d89edb9e548da48eaf09578f5b3a3aa215b18e3933141ff643b9b0e4865d27f
        var rhs = Math.pow(A*k_l_ratio**alpha-delta*k_l_ratio,-sigma)*A*Math.pow(k_l_ratio,alpha)*(1-alpha)/phi
        var log_rhs = Math.log(rhs)

        // objective = function(x) { return Math.pow(rhs - Math.pow((1-x),-eta)*Math.pow(x,sigma),2); }
        objective = function(x) { return Math.pow(log_rhs + eta*Math.log(1-x) - sigma*Math.log(x),2); }

        L0 = (sigma + Math.log(rhs))/(sigma+eta)

        var result = numeric.uncmin(objective,[0.1])

        var L = result['solution'][0]
        console.log("steady state L: ",L)
        console.log("objective at L: ",objective(L))
        console.log("objective at L: ",objective(0.39813726765942364))

        // var L = 1 - phi/((1-alpha)*A*Math.pow(k_l_ratio,alpha))
        var K = k_l_ratio*L
        var Y = A*Math.pow(K,alpha)*Math.pow(L,(1-alpha))
        var I = delta*K
        var C = Y - I
        var W = (1-alpha)*Y/L
        var R = alpha*Y/K
        
        // console.log("Steady state:")
        // console.log("A: ",A)
        // console.log("K: ",K)
        // console.log("C: ",C)
        // console.log("L: ",L)
        // console.log("I: ",I)
        // console.log("Y: ",Y)
        // console.log("W: ",W)
        // console.log("R: ",R)
        // console.log("K/L: ",k_l_ratio)

        // Construct coefficients for the log-linearized model
        var phi_1 = -(alpha+(eta-alpha)*L)/(1-L)
        var phi_2 = (alpha-1)*(1-beta*(1-delta))
        var phi_3 = (1-alpha)*(1-beta*(1-delta))
        var phi_4 = -rhoa*(1-beta*(1-delta))
        var phi_5 = A*Math.pow(K,alpha-1)*Math.pow(L,1-alpha)
        var phi_6 = alpha*A*Math.pow(K,alpha-1)*Math.pow(L,1-alpha) + 1 - delta
        var phi_7 = -C/K
        var phi_8 = (1-alpha)*A*Math.pow(K,alpha-1)*Math.pow(L,1-alpha)
        var phi_9 = phi_6-phi_8*alpha/phi_1
        var phi_10 = phi_7+phi_8*sigma/phi_1
        var phi_11 = phi_10*sigma/phi_1
        var phi_12 = (phi_9*sigma-phi_10*alpha)/phi_1
        var phi_13 = -phi_9*alpha/phi_1

        
        // Compute pi_4, pi_2, and pi_6
        var a = phi_3*phi_11 - sigma*phi_10
        var b = phi_2*phi_10-sigma*phi_9+phi_3*phi_12+sigma
        var c = phi_2*phi_9+phi_3*phi_13
        var pi_4 = (-b+Math.sqrt(Math.pow(b,2)-4*a*c))/2/a

        var pi_6 = (sigma*pi_4-alpha)/phi_1
        var pi_2 = phi_9 + phi_10*pi_4

        // Compute some additional coefficients to make things easier
        var phi_14 = phi_5 - phi_8/phi_1
        var phi_15 = phi_7 +  phi_8*sigma/phi_1
        var phi_16 = phi_2 - sigma*pi_4+phi_3*pi_6

        // document.write("Coefficients of linearized model:","<br><br>")
        // document.write("phi_1: ",phi_1,"<br>")
        // document.write("phi_2: ",phi_2,"<br>")
        // document.write("phi_3: ",phi_3,"<br>")
        // document.write("phi_4: ",phi_4,"<br>")
        // document.write("phi_5: ",phi_5,"<br>")
        // document.write("phi_6: ",phi_6,"<br>")
        // document.write("phi_7: ",phi_7,"<br>")
        // document.write("phi_8: ",phi_8,"<br>")
        // document.write("phi_9: ",phi_9,"<br>")
        // document.write("phi_10: ",phi_10,"<br>")
        // document.write("phi_11: ",phi_11,"<br>")
        // document.write("phi_12: ",phi_12,"<br>")
        // document.write("phi_13: ",phi_13,"<br>")
        // document.write("phi_14: ",phi_14,"<br>")
        // document.write("phi_15: ",phi_15,"<br>")
        // document.write("phi_16: ",phi_16,"<br>")
        // document.write("<br><br>")

        // Compute the reamining solution coefficients
        var numerator = phi_4 - phi_16*phi_14+phi_3*rhoa/phi_1
        var denominator = phi_16*phi_15 + sigma*(1-rhoa) + phi_3*rhoa*sigma/phi_1
        var pi_3 = numerator/denominator

        var pi_1 = phi_14 + phi_15*pi_3
        var pi_5 = (sigma*pi_3-1)/phi_1
        var pi_7 = Y/I*(1 + (1-alpha)*pi_5) - C/I*pi_3
        var pi_8 = Y/I*(alpha + (1-alpha)*pi_6) - C/I*pi_4
        var pi_9 = 1 + (1-alpha)*pi_5
        var pi_10 = alpha + (1-alpha)*pi_6

        var pi_11 = pi_9 - pi_5
        var pi_12 = pi_10 - pi_6

        var pi_13 = pi_9
        var pi_14 = pi_10 - 1

        // console.log("Solution coefficients:")
        // console.log("pi_1: ",pi_1)
        // console.log("pi_2: ",pi_2)
        // console.log("pi_3: ",pi_3)
        // console.log("pi_4: ",pi_4)
        // console.log("pi_5: ",pi_5)
        // console.log("pi_6: ",pi_6)
        // console.log("pi_7: ",pi_7)
        // console.log("pi_8: ",pi_8)
        // console.log("pi_9: ",pi_9)
        // console.log("pi_10: ",pi_10)
        // console.log("pi_11: ",pi_11)
        // console.log("pi_12: ",pi_12)
        // console.log("pi_13: ",pi_13)
        // console.log("pi_14: ",pi_14)

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
        var w_series = [0];
        var r_series = [0];

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
                    w_series.push(pi_11*a+pi_12*k);
                    r_series.push(pi_13*a+pi_14*k);

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
                w_series.push(pi_11*a+pi_12*k);
                r_series.push(pi_13*a+pi_14*k);

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
        sessionStorage.setItem('w_series',JSON.stringify(w_series));
        sessionStorage.setItem('r_series',JSON.stringify(r_series));
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
        $('#wProcess').highcharts({
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
                text: 'Wage',
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
                name: 'w',
                // color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*w_series[i])/100000,
                        })
                    }
                    return data;
                })()
            }
            ]
        });
        $('#rProcess').highcharts({
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
                text: 'Rental rate',
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
                name: 'r',
                // color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*r_series[i])/100000,
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
        var w_data = JSON.parse(sessionStorage.getItem('w_series'));
        var r_data = JSON.parse(sessionStorage.getItem('r_series'));
        var periods = parseInt(sessionStorage.getItem('periods'));

        var row1 = [];
        row1.push("period");
        row1.push("a");
        row1.push("k");
        row1.push("c");
        row1.push("y");
        row1.push("l");
        row1.push("i");
        row1.push("w");
        row1.push("r");

        var rows = [row1]
        for (i = 0; i <= periods; i++) {
            rows.push([i,a_data[i],k_data[i],c_data[i],y_data[i],l_data[i],i_data[i],w_data[i],r_data[i]]);
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