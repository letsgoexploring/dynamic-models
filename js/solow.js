$(document).ready(function () {


    $('#ParameterInput').on('submit', function (e) {
        e.preventDefault();
        var tmax = parseInt($('#tmax').val());
        var a1 = parseFloat($('#A1').val());
        var alpha1 = parseFloat($('#alpha1').val());
        var s1 = parseFloat($('#s1').val());
        var delta1 = parseFloat($('#delta1').val());
        var g1 = parseFloat($('#g1').val());
        var n1 = parseFloat($('#n1').val());
        var a2 = parseFloat($('#A2').val());
        var alpha2 = parseFloat($('#alpha2').val());
        var s2 = parseFloat($('#s2').val());
        var delta2 = parseFloat($('#delta2').val());
        var g2 = parseFloat($('#g2').val());
        var n2 = parseFloat($('#n2').val());                    
        var k0 = parseFloat($('#k0').val());
        var e0 = parseFloat($('#e0').val());
        var l0 = parseFloat($('#l0').val());
        var kstar1 = Math.pow(s1*a1/(delta1+n1+g1),1/(1-alpha1));
        var kstar2 = Math.pow(s2*a2/(delta2+n2+g2),1/(1-alpha2));
        var interval = 0;
        var enableMarks = true;
        var showCounterFactual = document.getElementById ("counterFactual");

        var input = document.getElementById ("steady");
            if (input.checked == true) {
                var k0 = kstar1;
            } else {
                var k0 = k0;
            };

        if (tmax<=30) {
            interval = 5;
        }
        else if (30<tmax&&tmax<=100) {
            interval = 20;
            enableMarks = false;
        }
        else if (100<tmax&&tmax<=500) {
            interval = 100;
            enableMarks = false;
        }
        else {
            interval = 500;
            enableMarks = false;
        };

        var y0 = A1*Math.pow(k0,alpha1);
        var c0 = (1-s1)*y0;
        var i0 = s1*y0;
        var K0 = e0*k0;
        var Y0 = e0*y0;
        var C0 = e0*c0;
        var I0 = e0*i0;

        var kSeries = [k0];
        var ySeries = [];
        var cSeries = [];
        var iSeries = [];
        var KSeries = [];
        var YSeries = [];
        var CSeries = [];
        var ISeries = [];
        var eSeries = [1];

        var kSeries2 = [k0];
        var ySeries2 = [];
        var cSeries2 = [];
        var iSeries2 = [];
        var KSeries2 = [];
        var YSeries2 = [];
        var CSeries2 = [];
        var ISeries2 = [];
        var eSeries2 = [1];
        // nSeries = [0];

        var k =k0;
        var kLag = k0;
        var eff=0;
        var effLag=0;
        var y=0;
        var k2 =k0;
        var kLag2 = k0;
        var eff2=0;
        var effLag2=0;

        for (i = 0; i <= tmax; i++) {
            if (i<5) {
                effLag = eSeries[eSeries.length-1];
                eff = (1+g1)*effLag;
                kLag = kSeries[kSeries.length-1];
                k = s1*a1*Math.pow(kLag,alpha1)+(1-delta1-n1-g1)*kLag;
                eSeries.push(eff);
                kSeries.push(k);
                ySeries.push(a1*Math.pow(kLag,alpha1));
                cSeries.push((1-s1)*a1*Math.pow(kLag,alpha1));
                iSeries.push(s1*a1*Math.pow(kLag,alpha1));
                KSeries.push(effLag*kLag);
                YSeries.push(effLag*a1*Math.pow(kLag,alpha1));
                CSeries.push(effLag*(1-s1)*a1*Math.pow(kLag,alpha1));
                ISeries.push(effLag*s1*a1*Math.pow(kLag,alpha1));

                effLag2 = eSeries2[eSeries2.length-1];
                eff2 = (1+g1)*effLag2;
                kLag2 = kSeries2[kSeries2.length-1];
                k2 = s1*a1*Math.pow(kLag2,alpha1)+(1-delta1-n1-g1)*kLag2;
                eSeries2.push(eff2);
                kSeries2.push(k2);
                ySeries2.push(a1*Math.pow(kLag2,alpha1));
                cSeries2.push((1-s1)*a1*Math.pow(kLag2,alpha1));
                iSeries2.push(s1*a1*Math.pow(kLag2,alpha1));
                KSeries2.push(effLag2*kLag2);
                YSeries2.push(effLag2*a1*Math.pow(kLag2,alpha1));
                CSeries2.push(effLag2*(1-s1)*a1*Math.pow(kLag2,alpha1));
                ISeries2.push(effLag2*s1*a1*Math.pow(kLag2,alpha1));
            }
            else {
                effLag = eSeries[eSeries.length-1];
                eff = (1+g2)*effLag;
                kLag = kSeries[kSeries.length-1];
                k = s2*a2*Math.pow(kLag,alpha2)+(1-delta2-n2-g2)*kLag;
                eSeries.push(eff);
                kSeries.push(k);
                ySeries.push(a2*Math.pow(kLag,alpha2));
                cSeries.push((1-s2)*a2*Math.pow(kLag,alpha2));
                iSeries.push(s2*a2*Math.pow(kLag,alpha2));
                KSeries.push(effLag*kLag);
                YSeries.push(effLag*a2*Math.pow(kLag,alpha2));
                CSeries.push(effLag*(1-s2)*a2*Math.pow(kLag,alpha2));
                ISeries.push(effLag*s2*a2*Math.pow(kLag,alpha2));

                effLag2 = eSeries2[eSeries2.length-1];
                eff2 = (1+g1)*effLag2;
                kLag2 = kSeries2[kSeries2.length-1];
                k2 = s1*a1*Math.pow(kLag2,alpha1)+(1-delta1-n1-g1)*kLag2;
                eSeries2.push(eff2);
                kSeries2.push(k2);
                ySeries2.push(a1*Math.pow(kLag2,alpha1));
                cSeries2.push((1-s1)*a1*Math.pow(kLag2,alpha1));
                iSeries2.push(s1*a1*Math.pow(kLag2,alpha1));
                KSeries2.push(effLag2*kLag2);
                YSeries2.push(effLag2*a1*Math.pow(kLag2,alpha1));
                CSeries2.push(effLag2*(1-s1)*a1*Math.pow(kLag2,alpha1));
                ISeries2.push(effLag2*s1*a1*Math.pow(kLag2,alpha1));


            }

        }

        sessionStorage.setItem('eSeries',JSON.stringify(eSeries));
        sessionStorage.setItem('kSeries',JSON.stringify(kSeries));
        sessionStorage.setItem('ySeries',JSON.stringify(ySeries));
        sessionStorage.setItem('cSeries',JSON.stringify(cSeries));
        sessionStorage.setItem('iSeries',JSON.stringify(iSeries));
        sessionStorage.setItem('KSeries',JSON.stringify(KSeries));
        sessionStorage.setItem('YSeries',JSON.stringify(YSeries));
        sessionStorage.setItem('CSeries',JSON.stringify(CSeries));
        sessionStorage.setItem('ISeries',JSON.stringify(ISeries));
        sessionStorage.setItem('periods',tmax);

        $('#capital_per_eff_worker').highcharts({
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
                text: 'Capital Per Effective Worker',
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
                name: 'k',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= tmax; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*kSeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'k original',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 5; i <= tmax; i++) {

                        if (showCounterFactual.checked===true) {
                            data.push({
                                x: i,
                                y: Math.round(100000*kSeries2[i])/100000,
                            })
                        }
                    }
                    return data;
                })()
            },
            ]
        });

        $('#output_per_eff_worker').highcharts({
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
                text: 'Output Per Effective Worker',
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
                name: 'y',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= tmax; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*ySeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'y original',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 5; i <= tmax; i++) {

                        if (showCounterFactual.checked===true) {
                            data.push({
                                x: i,
                                y: Math.round(100000*ySeries2[i])/100000,
                            })
                        }
                    }
                    return data;
                })()
            },
            ]
        });

        $('#consumption_per_eff_worker').highcharts({
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
                text: 'Consumption Per Effective Worker',
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
                name: 'c',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= tmax; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*cSeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'c original',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 5; i <= tmax; i++) {

                        if (showCounterFactual.checked===true) {
                            data.push({
                                x: i,
                                y: Math.round(100000*cSeries2[i])/100000,
                            })
                        }
                    }
                    return data;
                })()
            },
            ]
        });

        $('#investment_per_eff_worker').highcharts({
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
                text: 'Investment Per Effective Worker',
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
                name: 'i',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= tmax; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*iSeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'i original',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 5; i <= tmax; i++) {

                        if (showCounterFactual.checked===true) {
                            data.push({
                                x: i,
                                y: Math.round(100000*iSeries2[i])/100000,
                            })
                        }
                    }
                    return data;
                })()
            },
            ]
        });

        $('#capital_per_worker').highcharts({
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
                text: 'Capital Per Worker',
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
                data: (function () {
                    var data = [];

                    for (i = 0; i <= tmax; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*KSeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'k original',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 5; i <= tmax; i++) {

                        if (showCounterFactual.checked===true) {
                            data.push({
                                x: i,
                                y: Math.round(100000*KSeries2[i])/100000,
                            })
                        }
                    }
                    return data;
                })()
            },
            ]
        });

        $('#output_per_worker').highcharts({
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
                text: 'Output Per Worker',
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
                data: (function () {
                    var data = [];

                    for (i = 0; i <= tmax; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*YSeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'y original',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 5; i <= tmax; i++) {

                        if (showCounterFactual.checked===true) {
                            data.push({
                                x: i,
                                y: Math.round(100000*YSeries2[i])/100000,
                            })
                        }
                    }
                    return data;
                })()
            },
            ]
        });

        $('#consumption_per_worker').highcharts({
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
                text: 'Consumption Per Worker',
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
                data: (function () {
                    var data = [];

                    for (i = 0; i <= tmax; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*CSeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'c original',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 5; i <= tmax; i++) {

                        if (showCounterFactual.checked===true) {
                            data.push({
                                x: i,
                                y: Math.round(100000*CSeries2[i])/100000,
                            })
                        }
                    }
                    return data;
                })()
            },
            ]
        });

        $('#investment_per_worker').highcharts({
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
                text: 'Investment Per Worker',
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
                data: (function () {
                    var data = [];

                    for (i = 0; i <= tmax; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*ISeries[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'i original',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 5; i <= tmax; i++) {

                        if (showCounterFactual.checked===true) {
                            data.push({
                                x: i,
                                y: Math.round(100000*ISeries2[i])/100000,
                            })
                        }
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
        var k = JSON.parse(sessionStorage.getItem('kSeries'));
        var y = JSON.parse(sessionStorage.getItem('ySeries'));
        var c = JSON.parse(sessionStorage.getItem('cSeries'));
        var inv = JSON.parse(sessionStorage.getItem('iSeries'));
        var K = JSON.parse(sessionStorage.getItem('KSeries'));
        var Y = JSON.parse(sessionStorage.getItem('YSeries'));
        var C = JSON.parse(sessionStorage.getItem('CSeries'));
        var I = JSON.parse(sessionStorage.getItem('ISeries'));
        var periods = parseInt(sessionStorage.getItem('periods'));

        var row1 = [];
        row1.push("Period");
        row1.push("Capital per effective worker (K/EL)");
        row1.push("Output per effective worker (Y/EL)");
        row1.push("Consumption per effective worker (C/EL)");
        row1.push("Investment per effective worker (I/EL)");
        row1.push("");
        row1.push("Capital per worker (K/L)");
        row1.push("Output per worker (Y/L)");
        row1.push("Consumption per worker (C/L)");
        row1.push("Investment per worker (I/L)");
        row1.push("");
        row1.push("Labor efficiency (E)");

        var rows = [row1]
        for (i = 0; i <= periods; i++) {
            rows.push([i,k[i],y[i],c[i],inv[i],"",K[i],Y[i],C[i],I[i],"",e[i]]);
        }

        let csvContent = "data:text/csv;charset=utf-8,";
        rows.forEach(function(rowArray){
           let row = rowArray.join(",");
           csvContent += row + "\r\n"; // add carriage return
        }); 

        var encodedUri = encodeURI(csvContent);
        var link = document.createElement("a");
        link.setAttribute("href", encodedUri);
        link.setAttribute("download", "solow_data.csv");
        document.body.appendChild(link); // Required for FF

        link.click(); // This will download the data file named "solow_data.csv".
    }

}