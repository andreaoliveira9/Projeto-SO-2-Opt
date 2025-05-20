function plotTopology(Nodes,Links,selected)
% plotTopology(Nodes,Links,selected) - Plots the network topology with
%         the 'selected' nodes in red
%
% Nodes:    a matrix with 2 columns with the (x,y) coordinates of each node
% Links:    a matrix with 2 columns with the end nodes of each link
% selected: a row array with IDs of selected nodes (can be an empty row)

    nNodes= size(Nodes,1);
    %plot the links:
    plot([Nodes(Links(1,1),1) Nodes(Links(1,2),1)],[Nodes(Links(1,1),2) Nodes(Links(1,2),2)],'k-');
    hold on
    for i=2:size(Links,1)
        plot([Nodes(Links(i,1),1) Nodes(Links(i,2),1)],[Nodes(Links(i,1),2) Nodes(Links(i,2),2)],'k-')
    end
    clients= setdiff(1:nNodes,selected);
    %plot the non-selected nodes:
    plot(Nodes(clients,1),Nodes(clients,2),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10)
    for i=clients
        text(Nodes(i,1),Nodes(i,2),sprintf('%d',i),'HorizontalAlignment','center','Color','k','FontSize',6);
    end
    % plot the selected nodes:
    plot(Nodes(selected,1),Nodes(selected,2),'o','MarkerEdgeColor','r','MarkerFaceColor','w','MarkerSize',10)
    for i=selected
        text(Nodes(i,1),Nodes(i,2),sprintf('%d',i),'HorizontalAlignment','center','Color','r','FontSize',6);
    end    
    grid on
    hold off
end